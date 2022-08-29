# ISTRUZIONI

- _COMPILAZIONE:_

		make


- _ESECUZIONE:_

		./main.exe POT mu sigma walk

- _SIGNIFICATO E VALORI RAGIONEVOLI PER DATI INSERITI DA LINEA DI COMANDO:_
	- POT: 0 per potenziale armonico, 1 per il potenziale a doppia buca
	- mu: parametro della funzione d'onda di prova (ascissa (positiva e negativa) degli estremanti). Valori possibili su tutti i reali, sensati in (0,2).
  - sigma: parametro della funzione d'onda di prova (larghezza buca). Valori possibili su tutti i positivi, sensati in (0,0.8)
  - walk: passo nel metropolis. Valori possibili su tutti i positivi, sensati in (0.5,3) dalla regola 50/50.