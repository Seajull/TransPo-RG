``` python3
usage: tab1ToRef2.py [-h] [-f1 FASTA1] [-f2 FASTA2] [-ti TABINPUT] [-o OUT]
                     [-b FLANK] [-t TYPEF] [-i] [-v] [-w] [-c] [-te]

optional arguments:
  -h, --help            Affiche ces messages help
  -f1 FASTA1, --fasta1 FASTA1
                        Input de la référence fasta 1 (originelle).
  -f2 FASTA2, --fasta2 FASTA2
                        Input de la référence fasta 2 (nouvelle).
  -ti TABINPUT, --tabinput TABINPUT
                        Input du fichier tabulé associé à fasta 1.
  -o OUT, --output OUT  Fichier de sortie (format tabulé).
  -b FLANK              Taille de la fenêtre à sélectionner (par défaut : 50)
  -t TYPEF, --type TYPEF
                        Sélectionne uniquement les types souhaiter dans un
                        gff3 (exemple de format de l'option : 'mRNA,exon,cds'
                        (insensible a la casse)).
  -i, --index           Force l'indexation du fichier fasta 2.
  -v, --verbose         Active l'affichage
  -w, --warning         Désactive l'affichage des warnings.
  -c, --cds             Active la vérification des positions des CDS par
                        rapport aux gènes ou au mRNA.
  -te, --tempfile       Force la création du fichier fasta des séquences
                        sélectionnées, du ficher tabulé intermédiaire et du
                        fichier d'alignement .sam.
```
