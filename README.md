# LAMP Primer Auto design(LPAD) Tool

We developed the LAMP Primer Auto Design (LPAD) Tool to automate primer design for Loop-mediated Isothermal Amplification (LAMP). The tool scores primers based on key features such as melting temperature (Tm), GC content, terminal ΔG, primer distance, and secondary structure. Primer specificity is assessed using BLAST alignment against a user-defined reference genome. A Random Forest model is then used to predict primer quality and rank the primer sets. Evaluation on real and simulated datasets demonstrates a high model accuracy. In conclusion, LPAD provides a more efficient and comprehensive solution for LAMP primer design compared to existing tools.

## Table of Contents

- Environment
- Overall Pipeline
- Specific Implementation
- Ownership

## Environment
The environment configuration is provided in the file `environment.yml`.

### 1. **Operation System**
We strongly recommend using **Linux** or **MacOS** to run the tool.

### 2. **Dependencies**
  - python=3.8
  - primer3-py=2.0.3
  - biopython=1.78
  - pandas=2.0.3
  - blast=2.16.0
  - joblib=1.4.2
  - scipy=1.10.1
  - scikit-learn=1.3.0

## Overall Pipeline

<img src="https://github.com/user-attachments/assets/937c5124-fdd5-4e4f-b3fd-21793d9fd273" width="500" />

### 1. **Design Objective**：

LAMP amplification requires six primers: forward inner primer (FIP), backward inner primer (BIP), two outer primers (F3 and B3), and two loop primers (LF and LB). These primers are used to target specific sequence regions. Among them, loop primers are optional. This project mainly focuses on automatically designing the forward inner primer (FIP), backward inner primer (BIP), two outer primers (F3 and B3), and automatically outputting high-quality F1, F2, F3; B1, B2, B3 primers.

### 2. **Tool Workflow**：

- **Clone Respiratory**
    ```bash
    git clone https://github.com/YiCQi/BMI3-LPAD-Tool/
    cd BMI3-LPAD-Tool
    ```
    
- **Configure Environment**:
    - Ensure that Conda is installed on your computer. To install the necessary dependencies, please use the provided `environment.yml` file with Conda:

        ```bash
        conda env create -f environment.yml
        ```

    - This will create a new Conda environment with all the dependencies required for the project. If you want to use this program, please move to the current program's home directory and activate the conda environment:

        ```bash
        conda activate lpad
        ```

- **Input Data**：
    - The tool input includes alignment information of the target sequence (the region you want to amplify) and background genome information (the regions you do not want to amplify).

        ```bash
        # Use the default path (*** Note!!! This way requires the fasta file of the background genome 'hg38.fa' to be pre-placed in the folder './data/resource' ***)
        python LPAD.py
        ```

        ```bash
        # You can also specifies the custom path under the BMI3-LPAD-Tool folder
        python LPAD.py -i /path/to/your/input.fasta -r /path/to/your/reference.fasta -o /path/to/output/directory
        ```

    - If you want to specify the parameters for individual LAMP primers, modify the parameters for the corresponding primer file in the folder `'./configs'`

- **Sequence Alignment and Primer Generation**:
    - Primers are generated from the target sequence under certain constraints (length, spacing) and are strictly checked for specificity to prevent binding to the background (contaminating) genome
- **Defining Constraints and Scoring Criteria**:
    - Each candidate primer is evaluated against constraint conditions (such as Tm, GC content, free energy) and is assigned a score based on how well it meets these conditions.
    - Secondary structure prediction and dimer formation filtering are added, which further improves scoring accuracy and eliminates potential problematic primers.
- **Output Results**:
    - The tool will output a sorted list of primer combinations, ranked from highest to lowest score, providing the best candidates for user selection. The output files are stored in `'./data/output/Final_score/'`, and other intermediate files are stored in `'./data/output/Intermediate_file/'`

## Specific Implementation: 

### Feasibility of Algorithm: 

- **TM**:

Tm is estimated using the Nearest-Neighbor method. This method is currently considered to be the
 approximation method that gives the value closest to the actual value. 
 
The calculated Tm is affected by experimental conditions such as the salt concentration and oligo concentration,
 so it is preferred that Tm be calculated under fixed experimental conditions (oligo concentration at 0.1 µM, sodium
 ion concentration at 50 mM, magnesium ion concentration at 4 mM). 
 
The Tm for each region is designed to be about 65°C (64 - 66°C) for F1c and B1c, about 60°C (59 - 61°C) for F2,
 B2, F3, and B3, and about 65°C (64 - 66°C) for the loop primers. 
 
- **GC**%:
  
Primers are designed so that their GC content is between about 40% to 65%.

Primers with GC content between 50% and 60% tend to give relatively good primers. 

- **ΔG**:

The end of the primers serves as the starting point of the DNA synthesis and thus must have certain degree of
stability.   The 3’ ends of F2/B2, F3/B3, and LF/LB and the 5’ end of F1c/B1c are designed so that the free energy
is –4 kcal/ mol or less. The 5’ end of F1c after amplification corresponds to the 3’ end of F1, so that stability is 
important.

- **Primer Distance**:

<img src="https://github.com/user-attachments/assets/b17927cc-ba45-47dc-9bfc-23f14752519c" width="700" />

(Refer to the documentation of `Primer Explorer`)

- **Secondary Structure Prediction and Dimerization**:

Consider that there are no complementary sequences within the sequences (hairpin) and no complementary sequences between the sequences (dimer).

The design of hairpin: 1. The complementary segment length for hairpin structures should be between 6-12 bp. 2. The length of the loop segment should be within the range of 4-8 bp.

The design of dimer: Check whether there are 8-16bp complementary regions between the sequences.

## Ownership

Yihuan XU, Guobiao YE, Yicheng QI, Zhongyu SHI, and Hancheng YU

Zhejiang University-University of Edinburgh Institute, China.

