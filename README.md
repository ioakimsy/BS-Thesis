# Directory notes

## documents

- RM Slides
- SPP and Thesis templates
  
### SPP-2024

- Active LaTeX code for SPP 2024
  
### SPP-2024 drafts

- Drafts with comments for revision
  
### Thesis-Manuscript

- Active LaTeX code for BS thesis manuscript

## output

### 2D-Binary-PCA

- basic model:
  - Isotropic, homogeneous individual student learning rate, no memory
- Analysis
  - data
    - Contains a summary of characteristic variable ($m$) and time to learn ($T_{max}$ labeled as ttl) for different SA and class lengths (labelled as class_size)
      - $m$ and $T_{max}$ are in $\mu \pm \frac{\sigma}{\sqrt{N}}$ format
  - plots
    - Relevant plots for the summary of the data
    - $T_{max}$ vs $\lambda$ ($\lambda$ should be changed to $\rho$)
    - $m$ vs $\lambda$ ($\lambda$ should be changed to $\rho$)
    - $T_{max}$ vs $N$
- SA-$L$ 
  - $\rho_0$
    - trial #
      - animation
      - data
        - 2DBPCA-SA-$L$-$\rho_0$-data.csv
          - contains the state of each student at each generation
        - 2DBPCA-SA-$L$-$\rho_0$-fit_params.csv
          - contains fraction of learned per generation and power fit coefficients
          - power law: $y = p_1 \cdot x ^{p_2}$
      - images
        - 2DBPCA-SA-$L$-$\rho_0$-frame#.png
  
### 2D-Binary-PCA-IH

- modified basic model:
  - Isotropic, **inhomogeneous** individual student learning rate, no memory
- Analysis
  - data
  - plots

- SA-$L$
  - $\rho_0$-$\lambda_0$-$\delta\lambda$
    - trial #
      - animation (*not yet generated*)
      - data
        - 2DBPCA-SA-$L$-$\rho_0$-$\lambda_0$-$\delta\lambda$-data.csv
          - contains the state of each student at each generation
        - 2DBPCA-SA-$L$-$\rho_0$-$\lambda_0$-$\delta\lambda$-fit_params.csv
          - contains fraction of learned per generation and power fit coefficients
          - power law: $y = p_1 \cdot x ^{p_2}$
      - images (*not yet generated*)
        - 2DBPCA-SA-$L$-$\rho_0$-$\lambda_0$-$\delta\lambda$-frame#.png

### 1D models

- wolfram-ECAs
- wolfram-ECAs (circular boundary)
- PCA
  
## pluto-notebooks

- nothing important
  
## references

- references used for documents (SPP and BS thesis)
- bibtex not stored here
  
## scripts

- This is where the data comes from :0
  
## src

- Old and unused jupyter notebooked (.ipynb)