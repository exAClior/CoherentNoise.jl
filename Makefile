JL = julia


fmt:
	$(JL) -e 'using JuliaFormatter; format("src", style=BlueStyle())'