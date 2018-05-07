import argparse, sys, subprocess, tempfile


parser = argparse.ArgumentParser()

parser.add_argument("-g", "--getFlank",dest="getF", action="store_true", help="Lance uniquement la méthode getFlank() qui retourne leséquences flanquantes spécifiés dans le fichier tabulé et les stocks dans")
args = parser.parse_args()
parser.add_argument("-u", "--getuFlank",dest=str(args.getF)+"getu", action="store_true", help="Lance uniquement la méthode getFlank() qui retourne leséquences flanquantes spécifiés dans le fichier tabulé et les stocks dans")
args = parser.parse_args()
print(args.getF)
print(args)
