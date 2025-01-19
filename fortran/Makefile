all: interpol_modeles interpol_modeles_nlte

interpol_modeles: interpol_modeles.f
	gfortran -o interpol_modeles interpol_modeles.f

interpol_modeles_nlte: interpol_modeles_nlte_gfort.f
	gfortran -o interpol_modeles_nlte interpol_modeles_nlte_gfort.f

clean:
	rm interpol_modeles interpol_modeles_nlte
