Hello. Welcome to HIKE. In this repository, you will find necessary code, papers, etc. to get up and running. Look! We have a penguin too:  
https://na62.web.cern.ch/home/home.html


A) All the necessary papers, reviews and presentations in the 'Papaers_and_Reviews' folder (yes I'm aware of the typo, merely too lazy to correct it :) ).

B) 'dotroot_Data_Files' contains the data files containing the analysed reconstructions of various MonteCarlo simulations. Ask Rado for further details!

C) 'MC_Files' contains another README instructing how to run MC generation and Reconstuction using the NA62 framework. It also has a sample macro
    (K2pipi.mac) and the actual macro that has been in use since this project began (K02mumu.mac). Feel free to download them and make changes. You now
    control the fate of countless virtual Kaons, use this power wisely.

D) 'Codes' contains the bulk of what my time was spent doing and running over the course of the project. The .cpp files are faster on average and more
    robust in terms of handling errors, etc., so it is highly recommended to use .cpp files while running code using ROOT/C++.
    The 3 Jupyter Notebooks contain a startup guide to ROOT (ROOT Testing Ground), as well as the main notebooks (ROOTKaon and ROOTKaon_2), with
    instructions of their own within them. Note that everyhting has been coded in C++ and NOT PYTHON. I'm a C++ boi so I found a way to use C++ on
    Jupyter Notebooks, hence Python compilers will NOT BE ABLE to compile the code on there. Feel free to copy paste the code onto .cpp files for ease
    of use.


Please feel free to email me at (harrsh.goyal@students.iiserpune.ac.in) or (harrsh.goyal@cern.ch) or (scepxorus@gmail.com) for any issues or questions.
K0, enjoy.


P.S.
As an add-on, you can access the slides presentation (https://docs.google.com/presentation/d/1e2RKUplLKU9JjVrNqxO_XH1ZR0pCFdawO5WPVtZMP7o/edit?usp=sharing) containing a bunch of important results and basic information on how said results were obtained. Also, check the following google sheets doc (https://docs.google.com/spreadsheets/d/1hK2DYkHPSqPoH19BIojPw-5Nwnf47CUB4NFeC7NZwcs/edit?usp=sharing) which contains the data for various different mode of generation and fitting and how the fit and the uncertainity of the various parameters (CS, Cint and Phi0) turned out.
