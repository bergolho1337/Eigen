# Eigen

Examples programs that use the Eigen library.

# Pre-Requisites

  - Eigen 3.3.4. Available at: http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz
  - To install the Eigen library simple extract the file on a folder and on that folder do the following:
  ```sh
$ mkdir build; cd build
$ cmake ..
$ make
$ make install
```

> Remenber always to include the header (.h) on the folder where Eigen was installed by adding the -I on your Makefile (see the Example1 and Example2).

  - Some examples use the Qt4 library. To be able to run these you need to download the Qt4 libraries, which can be acquired by running on the terminal:

  ```sh
$ sudo apt-get install libqt4-dev
```
> For Qt4 programs you need to include the files needed by using the -I flag and also link the the libraries using the -L flag to your Makefile (see Example3).

 
