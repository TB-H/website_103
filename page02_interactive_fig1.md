---
layout: page
title: Predictive Regression 1
permalink: interactive_a
sidebar: true
interactive: predictive_regression_1.html
---
---

## Figure Description
Below is the predictive regression plot for both the exponential and linear models for bacterium 1 in which we have that
we have our samples along y and the linear to the left and the exponential to the right. 
We can notice how we appear to have greater overlap with the exponential cases when looking by eye at the predictive regression although we do have regions for both models in which data points are not within the 95% interval generated. 
We can especially note for the linear model how the discreptancies appear towards the ends of the line, noticing how for the right end as times increase, we begin to deviate from the linear model. We could note how the linear model could be seen as a first order Taylor approximation of the exponential model in which this would make sense since when doing the exponential model, we are fitting additional components of a series of polynomials. Let's see how this compares to the second bacterium.



<!-- The below line includes the interactive figure. Do not change! -->
<center>

{% include_relative interactives/{{page.interactive}} %}

</center>


