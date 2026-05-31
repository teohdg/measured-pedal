"""
analysis.py  -  Measured Pedal characterization toolkit
=======================================================
Computes gain, Bode frequency response, THD (via FFT), SNR, power/battery,
and measurement uncertainty from CSV data exported by your scope/analyzer.

It runs out-of-the-box: with no data files present it SYNTHESIZES realistic
sample data so you can verify your Python toolchain TODAY, before any hardware.

Usage:
    pip install numpy scipy matplotlib pandas
    python analysis.py            # runs the self-test demo, writes PNGs
Once you have real data, drop your CSVs next to this file and call the
functions directly (see the names at the bottom) or replace the demo paths.

CSV formats expected:
    waveform CSV : columns  time_s, voltage_v          (one scope capture)
    sweep CSV    : columns  freq_hz, Vin_rms, Vout_rms  (frequency response)
    drive CSV    : columns  drive_pct, Vin_rms, Vout_rms (gain vs knob)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.fft import rfft, rfftfreq


# Core math

def vpp_to_vrms(vpp):
    """Peak-to-peak -> RMS for a sine wave."""
    return vpp / (2.0 * np.sqrt(2.0))

def gain_db(vout_rms, vin_rms):
    """Voltage gain in dB."""
    return 20.0 * np.log10(np.asarray(vout_rms, float) / np.asarray(vin_rms, float))

def snr_db(vsig_rms, vnoise_rms):
    """Signal-to-noise ratio in dB."""
    return 20.0 * np.log10(vsig_rms / vnoise_rms)

def rms_from_waveform(t, v):
    """RMS of a captured waveform (DC removed)."""
    v = np.asarray(v, float) - np.mean(v)
    return np.sqrt(np.mean(v ** 2))

def thd_from_waveform(t, v, f0_hint=1000.0, n_harm=10):
    """
    Total Harmonic Distortion from a time-domain capture.
    THD = sqrt(V2^2 + ... + Vn^2) / V1, harmonic amplitudes read from the FFT.
    Returns (f0, thd_fraction, A1, harmonic_amps, freqs, spectrum).
    """
    t = np.asarray(t, float); v = np.asarray(v, float) - np.mean(v)
    fs = 1.0 / np.mean(np.diff(t))            # sample rate
    win = np.hanning(len(v))                  # window reduces spectral leakage
    spec = np.abs(rfft(v * win))
    freqs = rfftfreq(len(v), 1.0 / fs)

    band = (freqs > f0_hint * 0.5) & (freqs < f0_hint * 1.5)
    k1 = int(np.argmax(spec * band))
    f0 = freqs[k1]

    def amp_at(mult):
        k = int(np.argmin(np.abs(freqs - f0 * mult)))
        lo, hi = max(0, k - 2), min(len(spec), k + 3)
        return spec[lo:hi].max()

    A1 = amp_at(1)
    harms = np.array([amp_at(m) for m in range(2, n_harm + 1)])
    thd = float(np.sqrt(np.sum(harms ** 2)) / A1)
    return f0, thd, A1, harms, freqs, spec

def find_3db_corners(freq_hz, rel_db):
    """Return (low_corner, high_corner) where the relative response crosses -3 dB."""
    f = np.asarray(freq_hz, float); g = np.asarray(rel_db, float)
    def cross(order):
        idx = list(order)
        for a, b in zip(idx[:-1], idx[1:]):
            if (g[a] + 3) * (g[b] + 3) <= 0:
                tt = (-3 - g[a]) / (g[b] - g[a])
                return 10 ** (np.log10(f[a]) + tt * (np.log10(f[b]) - np.log10(f[a])))
        return None
    peak = int(np.argmax(g))
    return cross(range(0, peak + 1)), cross(range(len(f) - 1, peak - 1, -1))

def gain_db_with_uncertainty(vout, vin, rel_u_vout=0.02, rel_u_vin=0.02):
    """Gain in dB plus a Type-B uncertainty from instrument fractional errors (RSS)."""
    g = 20.0 * np.log10(vout / vin)
    rel_u = np.sqrt(rel_u_vout ** 2 + rel_u_vin ** 2)   # add in quadrature
    return g, 8.686 * rel_u                              # 8.686 = 20/ln(10)

def battery_life_hours(capacity_mAh, current_mA):
    return capacity_mAh / current_mA


# Plot helpers


def plot_bode(sweep_csv_or_df, ref_freq=1000.0, out="bode.png", conditions=""):
    d = sweep_csv_or_df if isinstance(sweep_csv_or_df, pd.DataFrame) else pd.read_csv(sweep_csv_or_df)
    d = d.copy()
    d["gain_dB"] = gain_db(d["Vout_rms"], d["Vin_rms"])
    i_ref = (d["freq_hz"] - ref_freq).abs().idxmin()
    d["rel_dB"] = d["gain_dB"] - d.loc[i_ref, "gain_dB"]
    low, high = find_3db_corners(d["freq_hz"].values, d["rel_dB"].values)
    plt.figure(figsize=(7, 4))
    plt.semilogx(d["freq_hz"], d["gain_dB"], "o-")
    if low:  plt.axvline(low,  ls="--", lw=1)
    if high: plt.axvline(high, ls="--", lw=1)
    plt.xlabel("Frequency (Hz)"); plt.ylabel("Gain (dB)")
    plt.title("Frequency Response  " + conditions)
    plt.grid(True, which="both", alpha=.4); plt.tight_layout(); plt.savefig(out, dpi=150)
    plt.close()
    return low, high

def plot_thd_curve(level_rms, thd_pct, out="thd_vs_level.png", conditions=""):
    plt.figure(figsize=(7, 4))
    plt.semilogx(level_rms, thd_pct, "o-")
    plt.xlabel("Input level (Vrms)"); plt.ylabel("THD (%)")
    plt.title("THD vs. Input Level  " + conditions)
    plt.grid(True, which="both", alpha=.4); plt.tight_layout(); plt.savefig(out, dpi=150)
    plt.close()

def plot_spectrum(freqs, spec, f0, thd, out="fft_spectrum.png"):
    plt.figure(figsize=(7, 4))
    plt.semilogy(freqs, spec / spec.max())
    plt.xlim(0, f0 * 10)
    plt.xlabel("Frequency (Hz)"); plt.ylabel("Normalized magnitude")
    plt.title(f"Output Spectrum (f0={f0:.0f} Hz, THD={thd*100:.1f}%)")
    plt.grid(True, which="both", alpha=.4); plt.tight_layout(); plt.savefig(out, dpi=150)
    plt.close()


# Synthetic data generators (so the toolchain is testable with NO hardware)


def synth_clipped_sine(f0=1000.0, fs=200000.0, dur=0.02, amp=2.0, clip=0.65, noise=0.001):
    """A soft-ish clipped sine to mimic a driven overdrive output."""
    t = np.arange(0, dur, 1.0 / fs)
    raw = amp * np.sin(2 * np.pi * f0 * t)
    v = np.clip(raw, -clip, clip)                  # hard clip stand-in for the demo
    v += noise * np.random.randn(len(t))           # add a little noise
    return t, v

def synth_sweep():
    """A Tube-Screamer-ish response: mid hump near ~720 Hz, HF rolloff."""
    f = np.array([20,30,50,80,100,200,300,500,800,1000,2000,3000,5000,8000,10000,15000,20000.])
    # crude model: low rolloff + mid emphasis + high rolloff
    g = (1/np.sqrt(1+(80/f)**2)) * (1+2.5/np.sqrt(1+((f-800)/600)**2)) * (1/np.sqrt(1+(f/8000)**2))
    vin = np.full_like(f, 0.05)
    vout = vin * g * 12
    return pd.DataFrame({"freq_hz": f, "Vin_rms": vin, "Vout_rms": vout})


# Demo / self-test


if __name__ == "__main__":
    print("Running self-test on synthetic data (no hardware needed)...\n")

    # 1) THD + spectrum from a synthetic driven waveform
    t, v = synth_clipped_sine(amp=2.0, clip=0.65)
    f0, thd, A1, harms, freqs, spec = thd_from_waveform(t, v)
    print(f"[THD]  f0 = {f0:.0f} Hz,  THD = {thd*100:.1f}%")
    plot_spectrum(freqs, spec, f0, thd)

    # 2) Bode from a synthetic sweep
    low, high = plot_bode(synth_sweep(), conditions="(demo data)")
    print(f"[Bode] -3 dB corners: low={low}, high={high} Hz")

    # 3) THD vs level curve (sweep clip threshold to emulate increasing drive)
    levels = np.array([0.01,0.02,0.05,0.1,0.2,0.5,1.0])
    thds = []
    for L in levels:
        tt, vv = synth_clipped_sine(amp=L*4, clip=0.65)
        _, th, *_ = thd_from_waveform(tt, vv)
        thds.append(th*100)
    plot_thd_curve(levels, thds, conditions="(demo data)")
    print(f"[THD curve] {[f'{x:.1f}%' for x in thds]}")

    # 4) Gain with uncertainty
    g, u = gain_db_with_uncertainty(2.0, 0.05)
    print(f"[Gain] {g:.1f} +/- {u:.2f} dB")

    # 5) SNR and battery life
    print(f"[SNR]  {snr_db(1.0, 0.30e-3):.1f} dB")
    print(f"[Batt] {battery_life_hours(500, 6):.0f} h at 6 mA on a 500 mAh cell")

    print("\nWrote: fft_spectrum.png, bode.png, thd_vs_level.png")
    print("Toolchain OK. Replace synthetic calls with your real CSVs when measuring.")
