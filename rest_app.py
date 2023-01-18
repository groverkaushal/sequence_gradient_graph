from flask import Flask, redirect, url_for, render_template, request
import requests
from requests.structures import CaseInsensitiveDict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logomaker as lm

app = Flask(__name__)

@app.route("/", methods=["POST","GET"])
def home():
    if request.method == "POST":
        inputseq = request.form["input_seq"]
        parameter_ = request.form["parameters"]
        colour = request.form["colours"]
        parameters_length = int(request.form["parameter_length"])
        options_select = request.form["options"]

        temp = ""
        parameters = {}
        for i in range(len(parameter_)):
            if parameter_[i] == " ":
                parameters[temp] = parameter_[i+1]
                temp = ""
                continue
            if parameter_[i-1] == " ":
                continue
            if parameter_[i] == "\n":
                continue
            if parameter_[i] == "\r":
                continue
            temp = temp + parameter_[i]
        temp = ""

        def fasta_to_list(sequence):
            sequence = sequence + ">\n"
            aligned_seq = []
            tt = ""

            for i in sequence:
                if i == "\r":
                    continue
                if i == ">":
                    subs = True
                if subs == True:
                    if i == "\n":
                        aligned_seq.append(tt)
                        tt = ""
                        subs = False
                        continue
                    else:
                        continue
                if subs == False:
                    if i == "\n":
                        continue
                    else:
                        tt = tt + i
            aligned_seq.pop(0)

            return aligned_seq

        seq = fasta_to_list(inputseq)

        colours = {"A":'Red',"T":'Blue',"G":'Green',"C":'Orange',"P":'Black',"Q":'Brown',"R":'Violet', "S":'Purple',"U":'Pink',"-":'White'}
        # colour = colour[0:len(colour)-1]
        colour = colour + "\r"
        temp = ""
        temp1 = ""
        for i in range(len(colour)-1):
            if colour[i] == " ":
                temp = colour[i-1]
                j = i 
                while colour[j+1] != "\r":
                    j = j + 1
                    temp1 = temp1 + colour[j]
                    
                colours[temp] = temp1
                temp = ""
                temp1 = ""

        # parameters = {"ATGCA":"P" , "TCATG":"Q" , "TAGTT":"R" , "GTACT":"S" , "GGAAT":"T"}

        def MSA_fasta(sequence):
            urladder = ">\n"

            url = "https://www.ebi.ac.uk/Tools/services/rest/muscle/run"

            headers = {
                'Accept': 'text/plain',
            }

            for i in range(len(sequence)):
                if i == len(sequence)-1:
                    urladder = urladder + sequence[i]
                else:
                    urladder = urladder + sequence[i] + "\n>\n"

            data = {
                'email': 'yehhaibot1@gmail.com',
                'format': 'fasta',
                'sequence': urladder,
            }

            resp = requests.post(url, headers=headers, data=data)

            a = resp.content.decode()

            url = "https://www.ebi.ac.uk/Tools/services/rest/muscle/status/" + a

            headers = CaseInsensitiveDict()
            headers["Accept"] = "text/plain"

            while resp.content.decode() != "FINISHED":
                resp = requests.get(url, headers=headers)
            
            url = "https://www.ebi.ac.uk/Tools/services/rest/muscle/result/" + a + "/aln-fasta"

            headers = CaseInsensitiveDict()
            headers["Accept"] = "text/plain"

            resp = requests.get(url, headers=headers)

            alignedseq = resp.content.decode()

            return alignedseq

        if options_select == "Yes":
            alignedseq = MSA_fasta(seq)
            aligned_seq = fasta_to_list(alignedseq)
        else:
            aligned_seq = []
            count = []
            for i in range(len(seq)):
                    count.append(len(seq[i]))
            for i in range(len(seq)):
                    temp = ""
                    temp = temp + seq[i]
                    if (max(count)-len(temp))%2 == 0:
                        temp = int((max(count)-len(temp))/2)*"-" + temp + int((max(count)-len(temp))/2)*"-"
                    else:
                        temp = ((max(count)-len(temp))//2)*"-" + temp + (((max(count)-len(temp))//2)+1)*"-"
                    aligned_seq.append(temp)


        def p_generator(sequence):
            l = ""
            if parameters_length % 2 != 0:
                for i in range((parameters_length//2),len(sequence)-(parameters_length//2)):
                    if sequence[i-(parameters_length//2):i+(parameters_length//2)+1] in parameters:
                        l = l + parameters[sequence[i-(parameters_length//2):i+(parameters_length//2)+1]]
                    else:
                        l = l + sequence[i]
            else:
                for i in range((parameters_length//2)-1,len(sequence)-(parameters_length//2)-1):
                    if sequence[i-(parameters_length//2)+1:i+(parameters_length//2)+1] in parameters:
                        l = l + parameters[sequence[i-(parameters_length//2)+1:i+(parameters_length//2)+1]]
                    else:
                        l = l + sequence[i]
            return l



        aligned_seq_p = []
        for i in range(len(aligned_seq)):
            aligned_seq_p.append(p_generator(aligned_seq[i]))


        counts_mat = lm.alignment_to_matrix(sequences=aligned_seq_p, to_type='information', characters_to_ignore='.-X')
        counts_mat.head()
        lm.Logo(counts_mat)

        plt.savefig('static/plot.png', format='png', dpi=500)
        

        return render_template("result.html", a = aligned_seq , b = aligned_seq_p , c = colours ,d = parameters)
        
    else:
        return render_template("index.html")




if __name__ == "__main__":
    app.run(debug=False,host="0.0.0.0")
