
## Install Anaconda on remote server
To be able to run the script on a remote server you will need the correct version of python. In our versions, we have used python 3.8.5. Installing directly on the OS is not adviced, so install using some virtual environment, for example Anaconda. Here, we go trhough how to install with Anaconda.

Install using wget:
 - Go to [Anaconda archive](https://repo.anaconda.com/archive/) and pick the latest version of Anaconda compatible with linux (Linux-x86_64.sh)
 - Log into your remote server using ssh.
 - In your root folder, type:
  ```
  $ wget https://repo.anaconda.com/archive/[Anaconda_version.sh]
  ```

 - Replace Anaconda_version.sh with the one you picked from the archive.
 - Next install Anaconda (executing in your root folder):
  ```
  $ bash [Anaconda_version.sh]
  ```

 - Accept the license and confirm installation location.
 - Next, you will need to set you conda environment variable. Do this by typing
  ```
  $ echo 'export PATH=~/anaconda3/bin:$PATH' >> ~/.bashrc
  $ source ~/.bashrc
  ```

You should now be all set to start using Anaconda. Check documentation for how to run different commands. 

Essentiallly, you want to create an environment for example

```
$ conda create --name ldpc_master_thesis python=3.8.5
$ source activate ldpc_master_thesis
```

Check that you are in the correct environment and start installing using conda install. You should now have the correct version of python when running scripts.



## Run jobs using screen

To run script remotely, use "screen" on linux.
To start a new screen type:

```
$ screen -S [screen_name]
```

You can pick your own "screen_name". 

Then, start your job inside the screen, and if you want to keep it running even after you have logged out from ssh press Ctrl+a+d.

To go into your screen again, type

```
$ screen -R [screen_name]
``` 

If you have forgotten your screen name:


```
$ screen -ls
```
