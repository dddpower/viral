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
7. Assuming well-fitted model or estimated latent space, set semantic change and grammaticality respectively as,
    $$
    \Delta \mathbf{z}[\tilde{x}_i] \stackrel{\text{def}}{=} \|\hat{\mathbf{z}} - \hat{\mathbf{z}}[\tilde{x}_i]\|_1
    $$  
    $$
    p(\tilde{x_i} \mid \mathbf{x}) \stackrel{\text{def}}{=} \hat{p}(\tilde{x_i} \mid \hat{\mathbf{z}_i})
    $$
    where $\hat{\mathbf{z}} \stackrel{\text{def}}{=} \frac{1}{N} \sum \hat{\mathbf{z}_i}$ and $\hat{p}(x_i \mid \hat{\mathbf{z_i}}) \stackrel{\text{def}}{=} \text{softmax}(\mathbf{W}\hat{\mathbf{z_i}} + \mathbf{b})$, for some learned model weights and bias $\mathbf{W}$, $\mathbf{b}$.
    
8. Denote priority score of mutation $\tilde{x_i}$ as sum of semantic change and grammaticality terms with coeffient $\beta$:
    $$
    a(\tilde{x_i};\mathbf{x}) \stackrel{\text{def}}{=} \Delta \mathbf{z}[\tilde{x}_i] + \beta p(\tilde{x_i} \mid \mathbf{x}) = \|\hat{\mathbf{z}} - \hat{\mathbf{z}}[\tilde{x}_i]\|_1 + \beta \hat{p}(\tilde{x_i} \mid \hat{\mathbf{z}_i})
    $$
9. Finally CSCS score is derived by applying rank function to each terms.
$$
a'(\tilde{x}_i; \mathbf{x}) = \text{rank}(\Delta z[\tilde{x}_i]) + \beta \cdot \text{rank}(p(\tilde{x}_i | \mathbf{x}))
$$