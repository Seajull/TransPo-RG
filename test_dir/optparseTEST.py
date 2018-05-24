import argparse, sys, subprocess,re, tempfile


parser = argparse.ArgumentParser()

parser.add_argument("-g", "--getFlank",dest="getF", action="store_true", help="Lance uniquement la méthode getFlank() qui retourne leséquences flanquantes spécifiés dans le fichier tabulé et les stocks dans")

parser.add_argument("-t", "--type",dest="typeF", default=None, help="Sélectionne uniquement les types souhaiter dans un gff3 (exemple de format de l'option : 'mRNA,exon,cds' (insensible a la casse)).")
args = parser.parse_args()
print(args.typeF)

print(re.findall("[a-zA-Z0-9]+",args.typeF))
typ=(re.findall("[a-zA-Z0-9]+",args.typeF))
o=[l.lower() for l in typ]
print(o)

