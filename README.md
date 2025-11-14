# Lipinski Check
A quick test for Lipinski's Rule of Five for all common molecular structure inputs and file types, including hyperlinks. 

## Table of Contents
* [Introduction](#introduction)
* [General Info](#general-info)
* [Technologies](#technologies)
* [Launch](#launch)
* [Example of Use](#example-of-use)

### Introduction

This is a quick executable, along with an example shown in the notebook, that takes most molecular structure files or formats (including SMILES & SDF/MOL as a file or rich text) as an input for any molecule and tests whether it passes or fails Lipinski's Rule of Five. 

### General Info

Lipinski's Rule of Five: https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five

### Technologies
* RDKit

### Launch
* Python

Sample SMILES: CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O

Sample Hyperlink: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/5288826/record/SDF?record_type=2d&response_type=display  

### Example of Use

Using ethanol as sample input in SMILES format:

Input: CCO 

Output:

<img width="720" height="322" alt="image" src="https://github.com/user-attachments/assets/11284dc3-0004-4768-b6bd-d84feb04ca63" />

