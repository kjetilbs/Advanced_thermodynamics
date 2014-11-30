#############################################################
# Iteration grid
#
# Using Julia to draw the iteration grid by outputting TikZ 
# commands to a .tex-file
#
# Author:     Kjetil Sonerud
# Updated:    2014-11-23 10:41:08
#############################################################

# Opening .tex-file
filename = "iterationGrid.tex"
f_handle = open(filename, "w")

# Defining macro for escaping string characters
macro R_str(s)
    s
end

# Writing preamble and start of document
#############################################################
preamble = "\\documentclass[tikz,border=5mm]{standalone}\n\\usepackage{tikz}\n\\usetikzlibrary{positioning}\n\n"
preamble = preamble * "\\begin{document}\n\\thispagestyle{empty}\n    \\begin{tikzpicture}\n"

write(f_handle, preamble)
#############################################################

# Writing the TikZ-picture
#############################################################
# Test
# tikzpicture = "        \\draw[thick,color=red] (0,1)--(1,2);\n"
# tikzpicture = tikzpicture * "        \\draw[thick,color=blue] (1,0)--(3,-1);\n"
# tikzpicture = tikzpicture * "        \\draw[thick,color=green] (0,0)--(1,1);\n"

# Creating indent
indent = "        "

tikzpicture = indent * R"\path [fill=white!90!gray] (0,0) rectangle (8.5,5.5);"*"\n"
tikzpicture = tikzpicture * indent * R"\draw [<->,>=stealth, ultra thick] (9,0) node[below]{$V \: \mathrm{[m^3]}$} -- (0,0) -- (0,6) node[left]{$T \: \mathrm{[K]}$};"*"\n"

# Initializing temperature and volume iteration array
tempIter    = [5:-1:1]
volumeIter  = [5:-1:1]

iterCounter = 0
addArrow    = 0.3

# Iterating on the volumes
for V in volumeIter
    # To achieve the desired temperature range; flip the 
    # range at every iteration
    tempIter = flipud(tempIter)
    addArrow = -addArrow

    # Iterating on the temperature
    for T in tempIter
        iterCounter += 1

        # Writing grid
        tikzpicture = tikzpicture * indent * R"\draw [fill=white!50!cyan] ("* string(V*1.5) * "," * string(T) * ") circle [radius=0.1] node at ("* string(V*1.5-0.4) * "," * string(T) * ") {"* string(iterCounter) * "};" * "\n"
    end
    tikzpicture = tikzpicture * indent * R"\draw [->,>=stealth, thick] ("* string(V*1.5-0.8) * "," * string(tempIter[1]+addArrow) * ") -- ("* string(V*1.5-0.8) * "," * string(tempIter[end]-addArrow) * ");" * "\n"
    tikzpicture = tikzpicture * indent * R"\draw [->,>=stealth, thick] ("* string(V*1.5-0.6) * "," * string(tempIter[end]-addArrow) * ") -- ("* string((V+1)*1.5-1) * "," * string(tempIter[end]-addArrow) * ");" * "\n"
end

write(f_handle, tikzpicture)
#############################################################

# Writing the postamble
#############################################################
postamble = "    \\end{tikzpicture}\n\\end{document}"

write(f_handle, postamble)
#############################################################

# Closing .tex-file
close(f_handle)

# Run pdflatex
run(`pdflatex iterationGrid.tex`)

# Clean up
run(`rm iterationGrid.log`)
run(`rm iterationGrid.aux`)



    # \begin{tikzpicture}
    #     \path [fill=white!85!gray] (0,0) rectangle (10,5);
    #     \draw [<->,>=stealth, ultra thick] (11,0) node[below]{$z \mathrm{[m]}$} -- (0,0) -- (0,6) node[left]{$r \mathrm{[m]}$};
    #     \foreach \y in {1, ...,5}
    #     {
    #         % Draw help lines 
    #         \draw [->, dashed] (0,\y) -- (10.5,\y);
            
    #         \foreach \x in {0, ..., 10}     
    #         {
    #             % Draw points
    #             \draw [fill=red] (\x, 0) circle [radius=0.1];
    #             \draw [fill=white!50!red] (\x, \y) circle [radius=0.1];
    #             \draw [fill=red] (\x, 5) circle [radius=0.1];
    #             %\draw[->] (1 - 0.1*\x,0.1*\y) -- (-0.1*\x + 0.2*\y,-0.1*\y);
    #         }
    #     }