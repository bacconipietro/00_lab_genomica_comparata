#!/bin/bash

OUTPUT="estrazione_risultati.tsv"

# Intestazione
echo -e "Pathway\tLikelihood_Data" > "$OUTPUT"

# Ciclo sulle cartelle K (da 1 a 5)
for k in {1..5}; do
    K_DIR="${k}K"
    
    # Definiamo il nome del file in base alla cartella K
    if [ "$k" -eq 1 ]; then
        FILENAME="Base_results.txt"
    else
        FILENAME="Gamma_results.txt"
    fi

    # Controllo esistenza cartella K
    if [ -d "$K_DIR" ]; then
        echo "Processando $K_DIR (cerco: $FILENAME)..."
        
        # Ciclo sulle cartelle N (da 1 a 10)
        for n in {1..10}; do
            TARGET_DIR="$K_DIR/${n}N"
            FULL_PATH="$TARGET_DIR/$FILENAME"
            
            if [ -f "$FULL_PATH" ]; then
                # Estrae la prima riga
                FIRST_LINE=$(head -n 1 "$FULL_PATH")
                # Scrive nel TSV
                echo -e "$TARGET_DIR\t$FIRST_LINE" >> "$OUTPUT"
            else
                # Opzionale: de-commenta la riga sotto se vuoi vedere gli errori
                # echo "File non trovato: $FULL_PATH"
                : # comando nullo (non fa niente)
            fi
        done
    else
        echo "Attenzione: Cartella $K_DIR non trovata."
    fi
done

echo "Finito! I dati sono stati salvati in $OUTPUT"
