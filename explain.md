# about data
used original cscs datatable as a source

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