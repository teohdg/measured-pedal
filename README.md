# Measured Pedal Project

Characterizing a soft-clipping analog gain stage, from prediction to silicon.

Tadeo Hernandez

I built the gain and clipping section of a Tube Screamer overdrive, but the pedal
was never really the point. The point was to run one small analog circuit through
the same loop an analog part goes through on its way to production: predict what it
should do on paper, simulate it, build it, measure it, and then see how badly the
real thing disagrees with the math.

The number I care about most: I measured 29.2% total harmonic distortion at 1 kHz.
The simulation predicted 28.77%. And the spectrum showed the odd harmonics beating
the even ones by about 36 to 1, which is what symmetric diode clipping is supposed
to do.

## What the circuit is

A non-inverting op-amp stage with two diodes wired backwards against each other in
the feedback loop. Quiet signals get amplified normally. Once the signal gets big
enough, the diodes start conducting and choke off the gain, so the peaks round over
instead of getting bigger. That rounding is the overdrive. A capacitor in the
gain-setting leg makes the amplification depend on frequency, which is where the
midrange emphasis comes from.

- TL072 op-amp, single 9 V supply with a mid-supply reference
- 1N4148 silicon clipping diodes in the feedback path
- Drive, tone, and level controls

## Results

| Parameter | Hand calc | LTspice | Measured |
|---|---|---|---|
| Peak gain (sim netlist) | 36.3 dB | 36.0 dB | -- |
| Pole frequency (sim netlist) | 720 Hz | 719 Hz | -- |
| Gain, minimum drive | 20.8 dB | -- | 23.4 dB |
| Gain, maximum drive | 41.4 dB | -- | 38.9 dB |
| THD at 1 kHz (clipping) | -- | 28.77% | 29.2% |
| Odd vs even harmonics | -- | predicted | 36 to 1 |

The [report](report/report.pdf) has the full data, the places where measurement and
prediction split apart, and the uncertainty math.

## How I worked

I calculated the gain, the corner frequencies, and the clipping threshold from the
component values first, before opening the simulator. Then LTspice, which matched
the hand math almost exactly: 36.0 dB against 36.3, and a 719 Hz pole against 720.
Then the breadboard, brought up one section at a time with a meter on every bias
point before I moved on. Then the bench measurements, then the comparison. The
comparison is the part that actually matters.

## What I found

The gain-versus-drive curve landed within about 2.6 dB of prediction at the ends of
the knob, but the measured range came out narrower than it should have. So I worked
backward from the numbers to figure out what feedback resistance they implied, and
the answer pointed at the potentiometer: its resistance never reached zero at one
end or its full value at the other. The disagreement was a worn control, not a
broken amplifier.

The frequency response came out the right shape, low in the bass, climbing, then
flat across the top. I measured it at minimum gain, because partway through testing
the drive pot quit and I swapped in a fixed resistor. The report works through why
the corner frequency shifted where it did, and why that shift lines up with the
minimum-gain condition rather than being a second, separate problem.

The distortion result is the one I am happiest with. Measured THD sat within half a
percent of the simulation, and the spectrum confirmed something I had only guessed
from the shape of the waveform, which is that symmetric clipping kills the even
harmonics and leaves the odd ones standing.

## The parts that went wrong

This is real bench work, so it includes the failures. I floated the supplies in the
first simulation and spent a while confused about why nothing amplified. On the
bench I measured the wrong node and read a gain 20 dB too low before I realized the
tone network was eating the signal. A probe scaling factor of ten fought me for an
afternoon. A resistor fell out. The drive pot developed an intermittent connection
and then died for good. The output oscillated until I traced it to a stability
capacitor I had left off the board. There is one input-level discrepancy I never
fully explained, so it is written up as open instead of dressed up as understood.
Anywhere I did not actually measure a value, it says TBD, not a guess.

## Layout

```
measured-pedal/
  README.md
  report/       report PDF and LaTeX source
  simulation/   LTspice files and output plots
  data/         measured CSVs
  analysis/     Python for gain, Bode, FFT and THD
  manual/       build guide and wiring reference
  photos/       breadboard and scope captures
```

## Tools

LTspice, Python with NumPy and Matplotlib, and a FNIRSI 2C53T that covers
oscilloscope, signal generator, and multimeter duty. Everything else was a
breadboard.

