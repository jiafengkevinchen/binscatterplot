# Binscatter implementation in Python

A Python wrapper of `binsreg` in R for binned scatterplots with automatic bandwidth selection and
nonparametric fitting (See [Cattaneo, Crump, Farrell, and Feng](https://arxiv.org/pdf/1902.09608.pdf)). 

- Uses `rpy2` and handles the input and output, so the user wouldn't have to worry about various R objects
- Uses Pythonic plotting capabilities like `seaborn` 
- Mimicks `seaborn.regplot` in usage 


## Minimal working example

![image](https://user-images.githubusercontent.com/24930289/65379164-aea61780-dc91-11e9-9d0b-f497d0d917ca.png)

# Author
Jiafeng Chen

