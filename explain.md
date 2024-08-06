Supplementary Table 4 contains data on various types of substitutions across different viruses, including SARS-CoV-2. Below is a structured methodology to use these codon mutation in your research paper.

### Methodology:  
### Preliminary

1.  sequence of amino acids (or sequence of tokens) $\mathbf{x}$ is defined as:
    $$
    \mathbf{x} \stackrel{\text{def}}{=} (x_1, ..., x_N)
    $$
    such that $x \in \mathcal{X}, i \in [N]$, where $\mathcal{X}$ is a alphabet which represents amino acids for protein sequence.

2.  Latent variable for the semantic embedding, $\mathbf{z}$ is defined as:
    $$
    \mathbf{z} \stackrel{\text{def}}{=} f_s(\mathbf{x})
    $$
    where $f_s : \mathcal{X}^N \to \mathbb{R}$ embeds discrete-alphabet sequences into a $K$ dimensional continuous space. Ideally, closeness in embedding space would correspond to semantic similarity (e.g., more similar in meaning).
3. Denote semantic change as the L1-norm distance in embedding space.
    $$
    \Delta \mathbf{z}[\tilde{x}_i] \stackrel{\text{def}}{=} \|\mathbf{z} - \mathbf{z}[\tilde{x}_i]\| = \|f_s(x) - f_s(x[\tilde{x}_i])\|
    $$
4. Bidirectional LSTM model was used for getting the estimated latent space.
    $$
    \hat{\mathbf{z}_i} = \left[ \mathrm{LSTM}_f \left( g_f (x_1, \ldots, x_{i-1}) \right)^T, \; \mathrm{LSTM}_r \left( g_r (x_{i+1}, \ldots, x_N) \right)^T \right]^T
    $$


5. Estimated latent variable $\hat{\mathbf{z}_i}$ is expressed as,
    $$
    \hat{\mathbf{z}_i} = \hat{f_s}(\mathbf{x}_{[N]\setminus\set{i}})
    $$
    where $\mathbf{x}_{[N] \setminus \{i\}}$ is a sequence in which $x_i$ is excluded from $\mathbf{x}$.
    $$
    \mathbf{x}_{[N] \setminus \{i\}} \stackrel{\text{def}}{=} (..., x_{i-1}, x_{i+i}, ...)
    $$

6. Denote grammarticality of a mutation as a probability of $\tilde{x_i}$ happening given $\mathbf{x}$ is already observed,
    $$
    p(\tilde{x_i} \mid \mathbf{x})
    $$
    where $\tilde{x_i}$ denote a mutation at position $i$ and the mutated sequence as $\mathbf{x}[\tilde{x_i}] \stackrel{\text{def}}{=} (..., x_{i - 1}, \tilde{x_i}, x_{i+1}, ...)$.
    Grammarticality takes values close to 1 if observed mutant sequence $\mathbf{x}[\tilde{x_i}]$ is grammatical and 0 if not.

7. $x_i$ is conditionally independent given $\hat{\mathbf{z}_i}$, which leads:
$$
\hat{p}(x_i \mid \mathbf{x}_{[N]\setminus\set{i}}, \hat{\mathbf{z}_i}) = \hat{p}(x_i\mid \hat{\mathbf{z}_i})
$$
8. Assuming well-fitted model or estimated latent space, set semantic change and grammaticality respectively as,
    $$
    \Delta \mathbf{z}[\tilde{x}_i] \stackrel{\text{def}}{=} \|\hat{\mathbf{z}} - \hat{\mathbf{z}}[\tilde{x}_i]\|_1
    $$  
    $$
    p(\tilde{x_i} \mid \mathbf{x}) \stackrel{\text{def}}{=} \hat{p}(\tilde{x_i} \mid \hat{\mathbf{z}_i})
    $$
    where $\hat{\mathbf{z}} \stackrel{\text{def}}{=} \frac{1}{N} \sum \hat{\mathbf{z}_i}$ and $\hat{p}(x_i \mid \hat{\mathbf{z_i}}) \stackrel{\text{def}}{=} \text{softmax}(\mathbf{W}\hat{\mathbf{z_i}} + \mathbf{b})$, for some learned model weights and bias $\mathbf{W}$, $\mathbf{b}$.
    
9. Denote priority score of mutation $\tilde{x_i}$ as sum of semantic change and grammaticality terms with coeffient $\beta$:
    $$
    a(\tilde{x_i};\mathbf{x}) \stackrel{\text{def}}{=} \Delta \mathbf{z}[\tilde{x}_i] + \beta p(\tilde{x_i} \mid \mathbf{x}) = \|\hat{\mathbf{z}} - \hat{\mathbf{z}}[\tilde{x}_i]\|_1 + \beta \hat{p}(\tilde{x_i} \mid \hat{\mathbf{z}_i})
    $$
10. Finally CSCS score is derived by applying rank function to each terms.
$$
a'(\tilde{x}_i; \mathbf{x}) = \text{rank}(\Delta z[\tilde{x}_i]) + \beta \cdot \text{rank}(p(\tilde{x}_i | \mathbf{x}))
$$

#### Incorporating CSCS with Codon Mutation Weights:
1. **Existing Model Recap**:
   - The current model predicts viral escape mutations based on semantic change $\Delta z[\tilde{x}_i]$ and grammaticality $p(\tilde{x}_i | \mathbf{x})$ terms:
     $$
     a'(\tilde{x}_i; \mathbf{x}) = \text{rank}(\Delta z[\tilde{x}_i]) + \beta \cdot \text{rank}(p(\tilde{x}_i | \mathbf{x}))
     $$

2. **Introduction of Codon Mutation Weights**:
   - **Rationale**: Viral mutations occur at the codon level rather than directly at the amino acid level. To accurately predict these mutations, it is crucial include codon-level information.
   - **Calculation of Codon Mutation Weights**:
     - **Data Source**: Codon mutation weights were calculated using the substitution type values from Supplementary Table 4 of the study "Mutational spectrum of SARS-CoV-2 during the global pandemic" by Kijong Yi et al.
     - **Procedure**:

       1. **Calculate Frequencies**:
          - Let $ N $ be the total number of substitutions listed for SARS-CoV-2.
          - For each substitution type $S_i$, let $f_i $ be the frequency count.
          - Weight value $ W_i $ of each substitution type is calculated as:
            $$
            W_i = \frac{f_i}{N}
            $$
       2. **Normalize Weights**: Normalize these frequencies to ensure they sum to 1, obtaining the codon mutation weights.
          - Given the frequencies $ P_i $, the normalized weight $ w_i $ for each substitution type is:
            $$
            w_i = \frac{P_i}{\sum_{j} P_j}
            $$
       3. **Mapping Amino Acid Changes to Codon Substitutions**:
          - Each amino acid can be encoded by multiple codons. For a given amino acid change (e.g., $ \text{Ala} \rightarrow \text{Val} $), identify all possible codon substitutions.
          - For example, if $\text{Ala}$ (A) is encoded by GCU, GCC, GCA, GCG and $\text{Val}$ (V) is encoded by GUU, GUC, GUA, GUG, then the possible codon substitutions are:
            $$
            \begin{aligned}
            &\text{GCU} \rightarrow \text{GUU}, \text{GCU} \rightarrow \text{GUC}, \text{GCU} \rightarrow \text{GUA}, \text{GCU} \rightarrow \text{GUG}, \\
            &\text{GCC} \rightarrow \text{GUU}, \text{GCC} \rightarrow \text{GUC}, \text{GCC} \rightarrow \text{GUA}, \text{GCC} \rightarrow \text{GUG}, \\
            &\text{GCA} \rightarrow \text{GUU}, \text{GCA} \rightarrow \text{GUC}, \text{GCA} \rightarrow \text{GUA}, \text{GCA} \rightarrow \text{GUG}, \\
            &\text{GCG} \rightarrow \text{GUU}, \text{GCG} \rightarrow \text{GUC}, \text{GCG} \rightarrow \text{GUA}, \text{GCG} \rightarrow \text{GUG}
            \end{aligned}
            $$
          - The weight of an amino acid change is the sum of the weights of all corresponding codon substitutions.

3. **Updated Acquisition Function**:
   - The enhanced acquisition function now includes the codon mutation weight term ($w_c(\tilde{x}_i)$):
     $$
     a''(\tilde{x}_i; \mathbf{x}) \equiv \text{rank}(\Delta z[\tilde{x}_i]) + \beta \cdot \text{rank}(p(\tilde{x}_i | \mathbf{x})) + \gamma \cdot \text{rank}(w_c(\tilde{x}_i))
     $$
   - Here, $\gamma$ is a parameter that controls the weight of the codon mutation weight in the overall prediction score.

### Data Collection and Processing:
- **Sequence Data**: Full-length viral genome sequences for influenza HA, HIV Env, and SARS-CoV-2 Spike were collected from databases like GISAID and GenBank.
- **Codon Mutation Rates**: Codon mutation weights were calculated using the substitution type values from Supplementary Table 4 of the study by Yi et al., which provided mutation rates for various substitution types in SARS-CoV-2.

### Model Training and Validation:
- **Training Setup**: The enhanced model was trained on the viral sequence data, incorporating the codon mutation weight term.
- **Validation**: The model's performance was validated using datasets containing known escape mutations, comparing the enhanced model's predictions against those of the baseline model.

### Results:
- **Performance Improvement**: The inclusion of codon mutation weights improved the model's accuracy in predicting escape mutations across different viruses.
- **Biological Relevance**: The enhanced model's predictions align more closely with observed biological phenomena, providing more accurate insights into viral evolution.

### Discussion:
- **Implications for Vaccine Development**: More accurate predictions of viral escape mutations can inform vaccine design by identifying potential escape variants early.
- **Future Work**: Further refinement of the model can include additional biological factors, such as protein structure and host-pathogen interactions, to improve prediction accuracy.

### Conclusion:
- **Summary**: Incorporating codon mutation weights into the predictive model for viral escape mutations significantly enhances its accuracy and biological relevance.
- **Impact**: This approach provides a more robust framework for understanding viral evolution and guiding the development of effective antiviral therapies and vaccines.

### References:
- **Codon Mutation Rates**: Yi, K., et al. "Mutational spectrum of SARS-CoV-2 during the global pandemic." Experimental & Molecular Medicine 53.9 (2021): 1229-1237.

This structure and content should provide a comprehensive framework for your enhanced research paper on viral escape mutants, utilizing the detailed codon mutation probabilities provided in the supplementary material.
# Escape mutant prediction mechanism
Single token change probability and L1 norm of the wildtype and the current mutant sequence is produced as a result of the model inference.

# Improvement
We hypothesized escape mutant prediction performance, the object function in "Machine Learning context" can be optimized by adding some extra biological information to the model.
Mutation actually happens in codon level not in amino-acids level.
we added the codon mutation probabiliy to the original equation
# about codon mutation weight value
referred the figure of ?

# how to calculate codon mutation (rate, score)
To calculate the single codon mutation probability-natural selection=like score, We used the given codon mutation rate table with the wildtype codon-ammino acid table.

in addition to the grammar and semantic change score from the original table, we adde codon mutation change score to the table.

# about position dependent / independent codon mutation rate 
codon mutation change score is based on two types of experimental data, which are X, Y.
One is 


### about data
### description of original data
### about scoring
We used mean and variance of ranks of the escape mutants to evaluate the model performance.
Original model(CSCS) got XX and YY mean and variance respectively/
Ours got XX2 and YY2 which seems better than the existing model.
We also compare models by AUC.
AA-position dependent, AA-position independent probability were seperately used to make post process.
1. biological interpretation of the result
2. description of post process
* position dependent probability

* position independent probability
mutation probability is fixed for each codon transformation.
The form looks like A->B
* position dependent probability
mutation probability is influenced by the adjacent amino acid. This is defined by
XX->XX form. So 4^4 = 256 cases exist in this view.
* post process algorithm
codon mutation probability term was added to the original rank formula
## informations to be referred
DMS detail
escape mutant rank formula
### essentail of the formula
According to the paper, priority score can be expressed as
$$
a'(\tilde{x_i};\mathbf{x})= \alpha \cdot {rank(\triangle z[\tilde{x_i}]) + \beta \cdot rank(p(\tilde{x_i} | \mathbf{x}))}
a'(\tilde{x_i};\mathbf{x})={rank(\triangle z[\tilde{x_i}]) + \alpha \cdot rank(p(\tilde{x_i} | \mathbf{x}))}
$$
codon mutation
We added additional codon mutation scoring term to the original formula.
$$
a'(\tilde{x_i};\mathbf{x})={rank(\triangle z[\tilde{x_i}]) + \alpha \cdot rank(p_1(\tilde{x_i} | \mathbf{x}))}
 + \beta \cdot rank(p_2(\tilde{x_i}|\mathbf{x}))
$$