# Nonparametric research on Alzheimer disease

This is the repository for the project of the _Nonparametric Statistics_ course at Politecnico di Milano. 

The aim of this research is to detect the main factors that contribute to the development of the disease, answering three main questions:

- Is there statistical difference between demented and nondemented patients?
- Can the dementia diagnosis be predicted using a non-specific medical test?
- Which type of patients are more likeable to become demented?

The first part of the study consists in making some preliminary inference about the difference in distribution between health and disease patients, and spotting some outliers that influence the distributions itself.
Then statistical difference between the two groups was assessed through permutational ANOVA tests.

In order to answer the second question a nonparametric logistic regression model was fitted using cubic splines, showing the (very good!) results in terms of prediction and comparing them with an vanilla logistic model. 

Finally a survival analysis lead to a better knowledge about the development of the disease in time, demonstrating the impact of gender on patients' survival probability.

For a more detailed insight of the study and the mathematical reference please refer to the [report](https://github.com/edoardopalli/Alzheimer-Detection/blob/main/Report.pdf).

### References
Computation were performed using 
>  R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

## Authors
Edoardo Palli - Politecnico di Milano  \
Elena Musiari - Politecnico di Milano \
Francesca di Filippo - Politecnico di Milano \
Erica Manfrin - Politecnico di Milano \
Master of Science in Mathematical Engineering students at Politecnico di Milano.
