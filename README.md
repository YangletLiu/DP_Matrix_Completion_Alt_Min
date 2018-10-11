# DP_Matrix_Completion_Alt_Min
##Abstract
Data privacy is fundamentally important in big data analysis, machine learning and distributed systems. To guarantee practical privacy, the notion of differential privacy is adopted to maximize the accuracy of queries in statistical database while minimizing the possibility of identifying records. In the matrix completion problem, recovery accuracy and data privacy are two contradicting aspects, since revealing more entries will increase the estimation accuracy but sacrifices privacy. Thus, privacy-preserving matrix completion is a challenging issue. In this paper, we propose a novel scheme for differentially private matrix completion with provably better accuracy bound. First, we present our differential private version of the practically efficient alternating minimization algorithm. Secondly, we present theoretical results for ensuring required privacy. Thirdly, on real world data sets, we compare the performances of our scheme with that of state-of-the-art schemes.



##Project Summary
+ We use Matlab to achieve a practical approach to complete matrix. 
+ For synthetic data, we compare our Differentially Private Matrix Completion via Alternating Minimization with Private Frank-Wolf algorithm, since both are designed for low-rank matrix completion, to show the advantage of our algorithm. We conduct experiment to recover a low-rank matrix with rank-$r$, from observed elements in the subset $\Omega$. 
+ For recovery error, we adopt the relative square metric, defined as $\text{RSE}=\frac{\|\widehat{\bm{M}}-\bm{M}\|_F}{\|\bm{M}\|_F}$ and $\text{RMSE}=\sqrt{\frac{1}{m}\sum\limits_{i=1}^m(y_i-\hat{y}_i)^2}$.
+ For running time, varying the matrix size and fixing other parameters, we measure CPU time in seconds.
+ For convergence speed, we measure the decreasing rate of the RES across the iterations by linearly fitting the measured RSEs (in log scale). We include those plots due to two reasons: 1) both algorithms are iterative; 2) the decreasing speed of the RSE provides explanations for the observed performance of the recovery error and the running time.



##Data
As we want to preserve privacy of every user, and the output for each user is $n$-dimensional, we expect the private recommendations to be accurate only when $m\gg n$. Due to this constraint, we conduct experiments on the following datasets: 1) Synthetic: We generate a random rank-one matrix $\bm{M}=\bm{U}\bm{V}^T$ with unit $\ell_\infty$-norm, $m=500K$, and $n = 400$, 2) Jester: This dataset contains $n = 100$ jokes, and $m \approx 73K$ users, 3) MovieLens10M (Top 400): We pick the $n = 400$ most rated movies from the Movielens10M dataset, resulting in $m \approx 70K$ users, 4) Netﬂix (Top 400): We pick the $n = 400$ most rated movies from the Netflix prize dataset, resulting in $m \approx 474K$ users, and 5) Yahoo! Music (Top 400): We pick the $n = 400$ most rated songs from the Yahoo! music dataset, resulting in $m \approx 995K$ users. We rescale the ratings to be from $0$ to $5$ for Jester and Yahoo! Music.



##Algorithm
\begin{algorithm}
    \caption{Differentailly Private Matrix Completion via Alternating Minimization: the Cloud Server}
    \label{alg2}
    \begin{algorithmic}[1]
        \REQUIRE Matrix dimension $m, n$, matrix rank $r$, total number of iterations $L\in \mathbb{N}$, probability $\frac{1}{p}$.
        \STATE Received one column of $P_{\Omega}(M)$ from all $n$ nodes, which is noted as $P_{\Omega_j}(\bm{M})$.
        \STATE $P_{\Omega}(M) = \cup_{j\in n} P_{\Omega_j}(\bm{M})$
        \STATE $\Omega = \cup_{j\in n} \Omega_j$
        \STATE Set top-$r$ left singular vectors of $\frac{1}{{p}}P_{\Omega_0}(\bm{M})$ to $\bm{U}_0$.
        \STATE {\bf Clipping step}: Set all elements of $\bm{U}^0$ that have magnitude greater than $\frac{2\mu\sqrt{r}}{\sqrt{n}}$ to zero and orthonormalize the columns of $\bm{U}^0$,
        \FOR{$\ell=0, \cdots, L-1$}
            \STATE Broadcast $\bm{U}^{\ell}$ to all $n$ nodes.
            \STATE Received $\bm{U}^{\ell+1}$ from nodes $j$.
            \STATE $\bm{U}^{\ell+1}=\frac{1}{n} \sum \bm{U}^{\ell+1}$
        \ENDFOR
    \end{algorithmic}
\end{algorithm}


\begin{algorithm}
    \caption{Differentailly Private Matrix Completion via Alternating Minimization at Node $j$, $j\in [n]$}
    \label{alg3}
    \begin{algorithmic}[1]
        \REQUIRE observed set $\Omega_j$, one column of values $P_{\Omega}(\bm{M})$: $P_{\Omega_j}(\bm{M})$, matrix dimension $m, n$, matrix rank $r$, total number of iterations $L\in \mathbb{N}$.
        \STATE Upload $P_\Omega_j(\bm{M})$.
        \STATE Partition $\Omega_j$ into $2L+1$ subset $\Omega_j^0,\, \cdots, \,\Omega_j^{2L}$ with each element of $\Omega_j$ belonging to one of the $\Omega_j^{\ell}$ with equal probability (sampling with replacement)
        \FOR{$\ell=0, \cdots, L-1$}
            \STATE Download $\bm{U}$
            \STATE $ \bm{V}_j^{\ell+1} \leftarrow \text{argmin}_{\bm{V}_j\in {\mathbb{R}^{r \times 1}}} {\| P_{\Omega_j^{\ell+1}}(\bm{U}^{\ell}\bm{V}_j)-\bm{M}_j\|}_F^2$, $j\in [n]$,
            \STATE $ \bm{U}^{\ell} \leftarrow \text{argmin}_{\bm{U}\in {\mathbb{R}^{m \times r}}} {\|P_{\Omega_j^{L+\ell+1}}(\bm{U} \bm{V}_j^{\ell})-\bm{M}_j\|}_F^2$, $j\in [n]$,
            \STATE $ \bm{U}^{\ell} \leftarrow \bm{U}^{\ell}+\bm{W}_j$, $\bm{W}_j\sim \mathcal{N}(\mathbf{0}, \sigma^2\bm{I})$, where $\mathcal{N}(\mathbf{0}, \sigma^2\bm{I})$ denotes a randomly Gaussian Noise.
            \STATE Upload $\bm{U}$ to server.
        \ENDFOR
        \STATE \Return{$X = \bm{U}^{\ell}\bm{V}^{\ell}$}  
    \end{algorithmic}
\end{algorithm}


\begin{algorithm}
    \caption{DP-Smoothed alternating least squares (DP-SAltLS): the Cloud Server}
        \label{alg4}
        \begin{algorithmic}[1]
        \REQUIRE  numbers of iterations $L\in \mathbb{N}$, target dimension $k$, coherence parameter $\mu$.
        \STATE Received one column of $P_{\Omega}(\bm{M})$ and $\Omega$ from all $n$ nodes, which is noted as $P_{\Omega_j}(\bm{M})$ and $\Omega_j$.
        \STATE $P_{\Omega}(M) = \cup_{j\in n} P_{\Omega_j}(\bm{M})$
        \STATE $\Omega = \cup_{j\in n} \Omega_j$
        \STATE $(\Omega_0,\Omega') \leftarrow \rm{Split}(\Omega,2),(\Omega_1,\cdots,\Omega_L)     \leftarrow \rm{SPLIT}(\Omega',L)$
        \STATE Set top-$r$ left singular vectors of $\frac{1}{{p}}P_{\Omega_0}(\bm{M})$ to $\bm{U}_0$.
        \FOR{$\ell=0, \cdots, L-1$}
            \STATE Broadcast $\bm{U}$ to all $n$ nodes.
            \STATE Received $\bm{U}^{\ell+1}_{j}$ from nodes 
            \STATE $\bm{U}^{\ell+1}=\frac{1}{n} \sum \bm{U}^{\ell+1}$
        \ENDFOR
        \end{algorithmic}
\end{algorithm}



\begin{algorithm}
    \caption{DP-Smoothed alternating least squares (DP-SAltLS)at Node $j$, $j \in [n]$}
        \label{alg5}
        \begin{algorithmic}[1]
        \REQUIRE observed set of indices $\Omega \in [m]\times[n]$ of an unknown matrix $\bm{M}\in\mathbb{R}^{m\times n}$ with entries $P_\Omega(\bm{M})$, numbers of iterations $L\in \mathbb{N}$, error parameter $\varepsilon > 0$, target dimension $k$, coherence parameter $\mu$.
        \STATE Upload $P_\Omega(\bm{M})_j$.
        \FOR{$\ell=0, \cdots, L-1$}
            \STATE Download $\bm{U}^{\ell}$
            \STATE $\widehat{\bm{V}}_j^{\ell+1} \leftarrow \rm{MedianLS}(P_{\Omega_j^{\ell+1}}(\bm{M}), \Omega_j^{\ell+1}, \bm{U}^{\ell}, L, k)$
            \STATE $\bm{V}_j^{\ell+1} \leftarrow \rm{SmoothQR}(\bm{V}_{\ell+1},\varepsilon,\mu)$


            \STATE $\bm{U}^{\ell+1} \leftarrow \rm{MedianLS}(P_{\Omega_j^{\ell+1}}(\bm{M}), \Omega_j^{\ell+1}, \widehat{\bm{V}}_j^{\ell+1}, L, k)+\bm{W}_j$, $\bm{W}_j\sim \mathcal{N}(\mathbf{0}, \sigma^2\bm{I})$, where $\mathcal{N}(\mathbf{0}, \sigma^2\bm{I})$ denotes a randomly Gaussian Noise.
            \STATE $\widehat{\bm{U}}^{\ell+1} \leftarrow \rm{SmoothQR}(\bm{U}_{\ell+1},\varepsilon,\mu)$

            \STATE Upload $\widehat{\bm{U}}^{\ell+1}$ to server.
        \ENDFOR
        \ENSURE Pair of matrices $(\bm{U}^{L-1},\bm{V}_j^{L})$
        \end{algorithmic}
\end{algorithm}