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