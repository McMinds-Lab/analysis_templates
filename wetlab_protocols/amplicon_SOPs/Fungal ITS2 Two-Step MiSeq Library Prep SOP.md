---
output:
  pdf_document: default
  html_document: default
---

**16S Two-Step MiSeq Library Prep**

**Materials:**

Milli-Q Water

AccuStart II PCR ToughMix (2x)

Fwd Primer (ITS7_MAf, 10 &micro;M)

Rev Primer (ITS4_MAr, 10 &micro;M)

Index Primers (10 &micro;M)

Sample DNA

**Notes & Protocol:**

1.  Combine primers with ToughMix in the following ratios, with enough for one reaction per sample, some negative controls, and pipetting error:

| ---                        | Single Reaction (&micro;L) | Final Conc. |
|----------------------------|------------------------|-------------|
| AccuStart II ToughMix (2X) | 6.25                   | 1X          |
| ITS7_MAf                   | 1.25                   | 0.2 &micro;M    |
| ITS4_MAr                   | 1.25                   | 0.2 &micro;M    |
| MilliQ Water               | 3.25                   | \-          |

2.  Aliquot master mix into PCR tubes (12 &micro;L for this 12.5 &micro;L reaction).
3.  Add 0.5 &micro;L of sample DNA to each experimental PCR tube. Use Milli-Q water for negative controls.
4.  Mix and spin down, and run the thermocycler protocol "1st step amplicon":

| Temperature (°C) | Time   | ---              |
|------------------|--------|------------------|
| 94               | 3 min | ---              |
| 94               | 45 s   | start repeat 35X |
| 50               | 60 s   | ---              |
| 72               | 90s    | end repeat       |
| 72               | 10 min | ---              |

5.  Run 10 uL of all samples and negative controls on a 1% low melting-point agarose gel with the 10 &micro;L well comb. Add loading dye directly to the PCR tubes and mix by pipetting, then space samples every OTHER well (multichannel pipet should work for this in the big gel rig).
6.  Perform second PCR, this time adding 12.5 &micro;L ToughMix and 9.5 &micro;L water directly to each well, then 1 &micro;L of each of the appropriate index primers directly to each well, and swirling therein the pipette tip with PCR product from the first amplification. Run the thermocycler program "Adaptor \-\> Indexed":

| Temperature (°C) | Time   | ---              |
|------------------|--------|------------------|
| 94               | 3 min | ---              |
| 94               | 45 s   | start repeat 12X |
| 63               | 60 s   | ---              |
| 72               | 90s    | end repeat       |
| 72               | 10 min | ---              |

8.  Run on gel to verify amplification & purity (the intermediate amplicon can be diluted 1:25 and run alongside these products as a negative control).
9.  Pool samples in equal volumes.
10. Clean up samples using standard bead protocol.
