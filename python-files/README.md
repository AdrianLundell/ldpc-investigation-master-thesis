# Description
The workflow is as follows:

  - Go to [compute_dmc](compute_dmc) to genereate a discrete memoryless channel that defines your channel. It is of course possible to define new channels. Just make sure that the output matches the existing channels.
  - Go to [density_evolution](density_evolution) to optimze the waterfall region through density evolution. Make sure to use the files that you have generated in [compute_dmc](compute_dmc).
  - Go to [mm_qc_pega](mm_qc_pega) to optimize the error floor by using the MM-QC-PEGA algorithm. Make sure to use the files that you have generated in [density_evolution](density_evolution).
  - The QC-code generated in the previous step can be used in [aff3ct](../aff3ct) to simulate the performance of the constructed code.
  - The [analysis](analysis) folder can be used to analyze various aspects of the constructed codes.
  - The [plots](plots) folder has been used to generate the figures for the article, but can also be used for analysis.
  

## VS Code
If you are running the code in vs code place 

´´´
    {
        "name": "Python: Current File",
        "type": "python",
        "request": "launch",
        "program": "${file}",
        "console": "integratedTerminal",
        "justMyCode": true,
        "cwd": "${fileDirname}"
    },
´´´

in your launch.json file in order to be able to run the scripts properly
