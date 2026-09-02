Pakiet Python do analizy wyników z pomiarów proteomicznych 18-plex TMT.
Generuje volcano ploty i diagramy Venna dla porównań MUT vs WT (samce i samice).

Jak zaintalować?

-m pip install ścieżka\do\proteomika


Jak uruchomić?


proteomika --input "Ścieżka\Proteomika\nazwa_pliku.csv" --output "Ścieżka\Proteomika\wyniki"


Z niestandardowymi progami:

proteomika --input "Ścieżka\Proteomika\nazwa_pliku.csv" --output "Ścieżka\Proteomika\wyniki" --fc 0.58 --pvalue 0.05


Z profilowaniem:

proteomika --input nazwa_pliku.csv --output wyniki --profile


Opisy argumentów


--input Plik CSV z proteomiki
--output Katalog wyjściowy
--fc Próg log2FC standard 1.5
--pvalue" Próg p-value 0.05
--profile Generuj raport profiler

Opis plików wynikowych
"volcano.html" Volcano plot
"significant.csv" Białka istotne dla każdego porównania 
"venn_M_vs_F_up.png" Venn plot: białka up-regulowane M vs F 
"venn_M_vs_F_down.png" Venn plot: białka down-regulowane M vs F 
"report.txt" Raport z wyników istotnych statystycznieS
"profile_report.txt" Raport profilera

Test

python -m pytest tests
