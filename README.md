## Installation instructions

- Install devtools: `install.packages("devtools")`
- Install statmod: `install.packages("statmod")`

## Formula for calculating Brier scores for censored data

Taken from page 67 from [this master thesis](https://www.duo.uio.no/handle/10852/69480):

<img width="595" alt="Screenshot 2020-09-27 at 19 08 25" src="https://user-images.githubusercontent.com/5699893/94421831-6f54b580-0186-11eb-8aa2-5346257db412.png">

The code for calculating the Brier score for one observation is given in the file [brier_score.R](brier_score.R), in the function `brier_score`.

Use the function `brier_scores` to calculate the Brier scores at n observed time points. This function requires estimated survival probabilities.

Similarly, use the function `brier_scores_fht` to calculate the Brier scores for an FHT model given parameter vectors `y0` and `mu`. Not sure if this is useful.
