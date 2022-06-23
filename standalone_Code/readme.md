This folder contains code files from scEpiSearch standalone tool. <br>
This is for demonstration purpose of the code only. <br>
In order to run the standalone tool Please download searchProject_new.tar.gz from http://reggen.iiitd.edu.in:1207/episearch/standalone/

If user wants to run the tool from command line. Please follow the following script :<br>
For finding projection of Query cells, use commandline_search.py script. <br>
Command to run script is : <br>
python3 commandline_search.py <br>
Following arguments user can edit within the script as mentioned below : <br>
####################### INPUT THE ARGUMENTS HERE######################################################################
<br>
	########## Query count file location ####<br>
	query_file = './sample_queries/human/GM/query_gm.txt'<br>
	########## Query peak file location #### <br>
	chr_file = './sample_queries/human/GM/chr_gm.bed' <br>
	########## No of top results for each query cell <br>
	top_study = 5 <br>
	########## No of Top clusters to search from <br>
	cluster = 5			     # No of Clusters <br>
	########## Active genes(Top 500) - 1 , Poised genes(Top 2000) - 2 <br>
	active_poised = 1		     <br>
	########## Query Type : Human - 1 , Mouse - 2 , Cross-Species - 3 <br>
	query_type = 1                      # 1 - Human , 2 - Mouse , 3 - Cross-Species <br>
	########## Annotated/Unannotated : Annotated - 1 , Unannotated - 2 <br>
	anno_unanno = 2		 <br>
	########### Accurate - 1 , Faster - 2 <br>
	accurate_faster = 2 <br>
################################################################################################################

For embedding of Query cells from different batches/species, use commandline_embedding.py script. <br>
Command to run script is : <br>
python3 commandline_embedding.py <br>
Following arguments user can edit within the script as mentioned below : <br>
#######################################################################################################
		########## Number of Dataset Queries <br>
		data =  2   #No of datasets <br>
		########## Number of Top results for each query <br>
		rps = 5 <br>
		########## Active(Top 500 Genes) - 1 , Poised(Top 1000 Genes) - 2 <br>
		active_poised = 2 <br>
		########## Top clusters to search from  <br>
		cls = 5 <br>
		########## Accurate - 1 , Faster -2 <br>
		var1 = 2   # accurate_faster <br>
		########## Annotated Reference only - 1 , Unannotated Reference - 2  <br>
		anno_unanno_emb = 2		#annotated_unannotated <br>
		########## Query Count files location <br>
		query_file = ['./queries_embedding/query_gm_GSE68103_human.csv','./queries_embedding/query_bcells_mouse.txt'] <br>
		########## Queries Peak files location (In same order as above) <br>
		chr_file = ['./queries_embedding/gm_human_GSE68103.bed','./queries_embedding/chr_bcells_mouse.bed'] <br>
		########## Species for each dataset (Human-1,Mouse-2) <br>
		vars = [1,2]		#species info <br>
############################################################################################################





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
