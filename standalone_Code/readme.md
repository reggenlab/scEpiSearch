This folder contains code files from scEpiSearch standalone tool. <br>
This is for demonstration purpose of the code only. <br>
In order to run the standalone tool Please download searchProject_new.tar.gz from http://reggen.iiitd.edu.in:1207/episearch/standalone/


User need to have python >=3.5+ install in their system. <br>

User can install all dependencies for the application using requirements.txt file present in the folder. <br>
One can use following command for same : <br>
python3 -m pip install --user -r /path_where_project/searchProject/requirements.txt <br>

or Following dependecies can be installed : <br>
MulticoreTSNE == 0.1 <br>
Pillow == 7.0.0 <br>
SimpSOM == 1.3.3 <br>
matplotlib == 3.0.3 <br>
networkx == 2.4 <br>
numpy == 1.16.4 <br>
numpy_indexed == 0.3.5 <br>
pandas == 0.24.2 <br>
reportlab == 3.3.0 <br>
rpy2 == 3.0.1 <br>
scikit_learn == 0.22.2.post1 <br>
scipy == 1.1.0 <br>
seaborn == 0.9.0 <br>
wordcloud == 1.7.0 <br>

User also need to have R installed in their system beforehand. <br>
R related packages required : GenomicRanges <br>
Rpy2 library used in this project creates a bridge between R and python. <br>

In order to install tkinter which provides GUI for the app, user needs to install it using following command : <br>
apt-get install python3-tk <br>

Finally user can run the application by reaching into the searchProject folder and running following command : <br>
"python3 scepisearch.py" <br>
If all dependencies work well, user can see the GUI. <br>

Note : If user runs into the error "OSError: cannot load library '/home/cell/R/lib/R/lib/libR.so': libBblas.so: shared Cannot open object file: there is no such file or directory " <br>
Solution : User needs to add path : export LD_LIBRARY_PATH = "/path-to-R-installation/lib/R/lib:$LD_LIBRARY_PATH"  <br>
------------------------------------------------------------------------------------------- ------------------------ <br>
If user is executing from location machine, GUI will work fine. <br>
Otherwise if user is executing from remote server, Following process has to be followed : <br>
1. If local system is windows, User needs to install Xming from : https://sourceforge.net/projects/xming/ <br>
	Afterwards, user should login to remote server from putty and tick X11 forwarding : <br>
	
2. If user is accessing remote server from linux based local system, then user can login to remote server with -X ssh login. <br>
For eg,  ssh -X username@192.168.X.X <br>
