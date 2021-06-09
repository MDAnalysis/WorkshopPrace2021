#!/bin/bash

for practical in $(ls ../*Practical/*ipynb)
do
  jupyter nbconvert --clear-output --inplace ${practical}
done
