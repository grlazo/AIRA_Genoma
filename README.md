## Genoma Genome Feature Finder (AIRA_Genoma)

The Genoma Genome Feature Finder (AIRA_Genoma) needs a few things to get started.

Basically you need the Ollama AI framework (****ollama.ai****) installed, 
and the anaconda (or miniconda3) to set up a python environment.

Ollama can be istalled on Linux, Windows, or Apple (check your requirements).

Once you have Ollama installed, you need a large language model (LLM). 
You can visit the website to discover all that is available. It's best to 
install your ollama instance under the conda environment to make sure both 
are working under the same environment.

For instance on a Windows machine you would start your Anaconda Powershell
and run the following commands:
```
$ ollama pull llama3.1:latest
$ ollama pull deepseek-r1:14b
$ ollama pull qwen2.5:7b
```
The above will get you started; you may change these later. Three LLMs 
pre-configured; you don't need all (edit  wheat_genome_browser.py if needed). 

You can check if ollama is running on your local machine by opening a web 
browser to the following address: http://localhost:11434/
It should display: ****ollama is running****
 
Next you will need to set up your conda environment; start by naming your
environment:
```
$ conda create --name Genoma python==3.12
$ conda activate Genoma
$ conda list
```
You will install needed packages in this environment; updates often occur 
so validated packages and versions will be included in a requirements.txt 
file. I like to install each individually to make sure everything loads 
properly.
```
$ pip install pandas
$ pip install matplotlib
$ pip install numpy
$ pip install tabulate
$ pip install langchain
$ pip install langchain_ollama
$ pip install langchain_experimental
$ pip install streamlit
```
There should be no issues under a Linux environment, but the peculiar steps 
suggested were encountered when walking someone through this process for a 
Windows machine (not sure about Apple). It seems to be a permission issue between user/administrator accounts (software installs require administrator permissions).

Place your GFF3 files in a directory called '***data***'. Start with one, or a 
few, to start until you're comfortable with the limitations on your machine settings.

Once the PDFs are in place the directory structure should look like:
```
./AIRA_Genoma/data/T_monococcum_TA299_HC_1A.gff3
./AIRA_Genoma/wheat_genome_browser.py
./AIRA_Genoma/analysis_functions.py
./AIRA_Genoma/custom_code_editor.py
./AIRA_Genoma/README.md_
./AIRA_Genoma/requirements.txt
```

Use streamlit to view the browser interface, start it with:
```
$ streamlit run wheat_genome_browser.py
```
Then open a browser with the link provided (usually http://localhost:5301/). 

Note: If you're in a Linux command shell, do:
```
$ grep '\$' README.md to see commands to issue after installing ollama.
```
