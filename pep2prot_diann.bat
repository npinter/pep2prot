@echo off
title pep2prot - DIA-NN edition
SET pep2prot_path=%~dp0
call C:\ProgramData\Miniconda3\Scripts\activate.bat D:\data\conda\pep2prot
python %pep2prot_path%pep2prot_diann.py
pause