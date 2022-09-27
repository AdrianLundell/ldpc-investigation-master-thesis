# Error Correction in NAND Flash Memories Using Low-Density Parity-Check Codes

This project enables the design of low-density parity-check codes (LDPC codes) of quasi-cyclicstructure for error correction of a chennel representing a NAND flash memory, achieving a maximum decoded bit error rate around 2 · 10−8 for a 1% raw bit error rate for a (9, 73, 256) code of code rate 0.88 over the BI-AWGN channel. This is accomplished by optimizing the code ensemble using density evolution to achieve a good waterfall threshold and minimizing trapping sets with a progressive edge-growth algorithm to lower the error floor of the decoding performance curve. The software tools developed is available in three [Python scripts](python-files) for code design and and a modified version of the C++ library (aff3ct)[aff3ct] for simulation (forked from [here](https://github.com/aff3ct/aff3ct)).

Link to our article [here]().
