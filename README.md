# Python for biologists project 
## Analysing processing body numbers and sizes in spermatogonial stem cells
    The aim of this project is to eliminate the need for manual analysis of
    individual image stacks and or excel files one by one, it allows the user 
    to get quick statistical results and plots.
### Project background
    Spermatogonial stem cells are cells that retain their ability to differenciate 
    and this can be tracked by a unique marker called cKIT.
    
 [cKIT](https://en.wikipedia.org/wiki/KIT_(gene))
 
    Spermatogonial cells positive for cKIT are demed as destined to undergo differentiation whereas
    cKIT negative cells retain their stem cell ability within this heterogenous cell group. This short 
    project aims to describe the translational represion differences between these two cell types by 
    analysing their processing bodies. Processing bodies are RNP granules that are indicators of 
    translationally represed mRNAs that are "trapped" within the cytoplasm.  
    
[processing bodies](https://pubs.acs.org/doi/10.1021/acs.biochem.7b01162)

***find out more***
[article1](https://www.frontiersin.org/articles/10.3389/fendo.2022.895528/full)
[article2](https://www.mdpi.com/2073-4409/9/3/745)

# How to run
Make sure you have Git and Conda installed before following these steps

1. Clone the git repository
    git clone <repository_url>
2. Navigate the directory
    cd project-directory
3. Create a Conda environment(if it does not already exist)
    conda create -n myenv python=3.8
4. Activate Conda environment
    conda activate myenv
5. Install `requirements.txt` file with pip
    pip install -r requirements.txt
6. Run the code using Python
    python code.py
   
7.You should have an output.csv file as your final result

