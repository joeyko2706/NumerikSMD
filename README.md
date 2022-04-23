# NumerikSMD

In diesem Repo sind die Dateien zu den Kursen Numerik und SMD (bisher nur SMD-A) aus dem Sommersemester 2022.

In Basisordner liegt der LaTeX-Header, Literaturdatenbank (`lit.bib`) sowie die der LaTeX-Header für matplotlib
und die matplotlibrc.
Außerdem haben wir in `programme.bib` die korrekten Quellen für die verwendete Software angegeben.

In dem Unterordner `vXXX` liegt dann ein Template für einen einzelnen Versuch. Der muss natürlich für die einzelnen
Abgaben in SMD umbeschrieben werden ^^.

## Mit Git arbeiten
https://guides.github.com/introduction/flow/.

1. Create a new branch using `git branch <name>`
2. Switch to it using `git checkout <name>`
3. Make changes and commit
4. Push the Branch using git `push -u origin <name>`
5. Open a Pull Request on github.

- Falls eine zu große Datei hochgeladen wurde und somit ein Fehler ausgegeben wird, oder wenn einfach der letzte 
    lokale commit gelöscht werden soll: `# git reset --soft HEAD^`