``` python3
usage: tab1ToRef2.py [-f1 FASTA1] [-f2 FASTA2] [-ti TABINPUT] [-b FLANK] [-c]
                     [-h] [-i] [-o OUT] [-t TYPEF] [-te] [-v] [-w]

Arguments requis :
  -f1 FASTA1, --fasta1 FASTA1
                        Input de la référence fasta 1 (originelle).
  -f2 FASTA2, --fasta2 FASTA2
                        Input de la référence fasta 2 (nouvelle).
  -ti TABINPUT, --tabinput TABINPUT
                        Input du fichier tabulé associé à fasta1.

Arguments optionnels :
  -b FLANK              Taille de la fenêtre à sélectionner (par défaut : 50)
  -c, --cds             Active la vérification des positions des CDS par
                        rapport aux gènes ou au mRNA.
  -h, --help            Affiche ces messages help
  -i, --index           Force l'indexation du fichier fasta 2.
  -o OUT, --output OUT  Fichier de sortie (format tabulé).
  -t TYPEF, --type TYPEF
                        Sélectionne uniquement les types souhaiter dans un
                        gff3 (exemple de format de l'option : "mRNA exon cds"
                        (insensible a la casse)).
  -te, --tempfile       Force la création du fichier fasta des séquences
                        sélectionnées, du ficher tabulé intermédiaire et du
                        fichier d'alignement .sam.
  -v, --verbose         Active l'affichage
  -w, --warning         Désactive l'affichage des warnings.
```
usage: tab1ToRef2.py [-f1 FASTA1] [-f2 FASTA2] [-ti TABINPUT] [-b FLANK] [-c]
                     [-h] [-i] [-o OUT] [-t TYPEF] [-te] [-v] [-w]

Arguments requis :
  -f1 FASTA1, --fasta1 FASTA1
                        Input de la référence fasta 1 (originelle).
  -f2 FASTA2, --fasta2 FASTA2
                        Input de la référence fasta 2 (nouvelle).
  -ti TABINPUT, --tabinput TABINPUT
                        Input du fichier tabulé associé à fasta1.

Arguments optionnels :
  -b FLANK              Taille des régions flanquantes à extraire de par et
                        d'autre de l'annotation (par défaut : 50).
  -c, --cds             Active la vérification des positions des CDS par
                        rapport aux gènes ou au mRNA.
  -h, --help            Affiche ces messages help
  -i, --index           Force l'indexation du fichier fasta 2.
  -o OUT, --output OUT  Fichier de sortie (format tabulé).
  -t TYPEF, --type TYPEF
                        Sélectionne uniquement les types souhaiter dans un
                        gff3 (exemple de format de l'option : "mRNA exon cds"
                        (insensible a la casse)).
  -te, --tempfile       Force la création du fichier fasta des séquences
                        sélectionnées, du ficher tabulé intermédiaire et du
                        fichier d'alignement .sam.
  -v, --verbose         Active l'affichage
  -w, --warning         Désactive l'affichage des warnings.
