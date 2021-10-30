# # ██▄      ▄▄▄▄▄   ▄████  ██   █    ▀▄    ▄ ▄▄▄▄▄▄   ▄███▄   █▄▄▄▄
# # █  █    █     ▀▄ █▀   ▀ █ █  █      █  █ ▀   ▄▄▀   █▀   ▀  █  ▄▀
# # █   █ ▄  ▀▀▀▀▄   █▀▀    █▄▄█ █       ▀█   ▄▀▀   ▄▀ ██▄▄    █▀▀▌
# # █  █   ▀▄▄▄▄▀    █      █  █ ███▄    █    ▀▀▀▀▀▀   █▄   ▄▀ █  █
# # ███▀              █        █     ▀ ▄▀              ▀███▀     █
# # ▀      █                                 ▀
# #
# # DSFalyzer
# #   Software to import, clean, and analyze DSF data from Roche LightCycler by
# #   second derivative Tm calling and isothermal analysis.
# # Michael Ronzetti, 2021

source('functions.R')
filePath <- './data/his_firstderiv.txt'

df.raw <- importData(filePath = filePath, type = 'FirstDeriv')
