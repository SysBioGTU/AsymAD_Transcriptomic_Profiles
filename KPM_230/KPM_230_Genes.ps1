

Set-Location "C:\Users\User\Desktop\KPM_230"

$filename = "C:\Users\User\Desktop\KPM_230\Binary_list.txt"

$path = "C:\Users\User\Desktop\KPM_230\results"

# Echo the file name
echo $filename

# Run the KPM analysis for the specific file
& java -Xmx1g -jar KPM-4.0.jar -graphFile="BIOGRID-ALL-4.4.236.tab3.txt(1).sif" -L1=0 -strategy=INES -algo=GREEDY -matrix1="$filename" -K=8

