import matplotlib.pyplot as plt
import tkinter
from tkinter import *


# auteur: Gijsbert Keja
# datum:30-01-2020
# studentenummer:618317


def main():
    """
    Roept alle functies aan en return de varaiabelen
    """
    gen_id_dict = openen()
    filter_gen_id = zoeken(gen_id_dict)
    gff_lijst = lezen_gff()
    gff_lijst_gefilterd = gff_gen_id(filter_gen_id, gff_lijst)
    grafiek(gff_lijst_gefilterd)


def openen():
    """
    opent het fasta bestand en maakt er 2 lijsten van een van de header waarmee alleen het id wordt meegegeven en
    een lijst van de proteines waarna deze de 2 lijsten in een dictionary zet
    :return: gen_id_dict een dict met als key de header en als value de proteïne die erbij hoort
    """
    try:
        fasta_openen = open('TAIR10_pep_20101214.fa', 'r')  # open de fasta file en leest deze alleen
        header = []
        proteines = []
        sequentie = ''
        for x in fasta_openen:
            x = x.strip()
            if '>' in x:
                if sequentie != '':
                    proteines.append(sequentie)
                    sequentie = ''
                code = x[1:]
                code = code.split(' ')
                header.append(code[0])
            else:
                sequentie += x.strip()
        proteines.append(sequentie)
        gen_id_dict = dict(zip(header, proteines))
        return gen_id_dict  # Dict met de header als key en de value als proteïnen
    except (IOError, FileNotFoundError, NameError):
        print('Bestand kan niet gelezen worden of gevonden worden of naam van var is fout')


def zoeken(gen_id_dict):
    """
    deze functie verandert de variabele van genID_dict en filter deze door middel van een regex
    :param: gen_id_dict een dict met als key de header en als value de proteïne die erbij hoort
    :return: filter_gen_id de dict die was ingegeven maar dan gefilterd op de eiwitconsensus
    """
    try:
        filter_gen_id = gen_id_dict  # geef de oude dic een nieuwe naam
        vervang_lijst = []
        for x in filter_gen_id:
            value = filter_gen_id[x]
            concensus_filter = re.search('[LIVMFYC].[HY].D[LIVMFY]K..N[LIVMFYCT][LIVMFYCT][LIVMFYCT]', value)
            if concensus_filter is None:
                vervang_lijst.append(x)
        for x in vervang_lijst:
            if x in vervang_lijst:
                del filter_gen_id[x]
        return filter_gen_id  # de oude dic gefilterd met de concensus patroon
    except (TypeError, NameError):
        print('string verwacht of naam van var is fout')


def lezen_gff():
    """'
    maak van ban het gff bestand wat op tabs gesplit is een 2d lijst
    :return: gff_lijst een 2d lijst van de .gff file waarbij gesplit is op de tabs
    """
    try:
        gff_lijst = list()
        gff = open('TAIR10_GFF3_genes.gff', 'r')
        for x in gff:
            x = x.replace('\n', '')  # haalt alle \n weg en vervangt ze
            entry = x.split('\t')
            gff_lijst.append(entry)
        return gff_lijst  # 2d lijst van de gff file
    except (IOError, FileNotFoundError, NameError):
        print('Bestand kan niet gelezen worden of kan niet gevonden worden of naam van var is fout')


def gff_gen_id(filter_gen_id, gff_lijst):
    """
    De volledige gff_file loopt over de de dict heen en als de key van de dic voorkomt in de lijst wordt deze
    toegevoegd aan de lege gff lijst
    :param: filter_gen_id de dict die was ingegeven maar dan gefilterd op de eiwitconsensus, gff_lijst een 2d lijst van
     de .gff file waarbij gesplit is op de tabs
    :return: gff_lijst_gefilterd een 2d lijst van de gff file maar dan met de id's die overeekomen met de concensus
    patroon eiwitten
    """
    try:
        gff_lijst_gefilterd = []
        for i in gff_lijst:
            if i[-1] in filter_gen_id:  # als de laatste kolom gelijk is aan de key gaat die naar de volgende stap
                gff_lijst_gefilterd.append(i)
        return gff_lijst_gefilterd  # de 2d gff file is gefilterd
    except (NameError, KeyError):
        print('Je key komt niet voor in je dict of naam van var is fout')


def grafiek(gff_lijst_gefilterd):
    """'
    maak een grafiek van de 2d lijst van de gff file hij telt de aantal gene per chromosoom op waarna deze
    in een lijt wordt gezet voor de grafiek
    :param:gff_lijst_gefilterd een 2d lijst van de gff file maar dan met de id's die overeekomen met de concensus
    patroon eiwitten
    """
    try:
        chr1 = 0
        chr2 = 0
        chr3 = 0
        chr4 = 0
        chr5 = 0
        chromosomen_nummers = ['1', '2', '3', '4', '5']  # geven de staven het juiste nummer
        for x in gff_lijst_gefilterd:
            if x[0] == 'Chr1':
                chr5 += 1
            if x[0] == 'Chr2':
                chr5 += 1
            if x[0] == 'Chr3':
                chr5 += 1
            if x[0] == 'Chr4':
                chr5 += 1
            if x[0] == 'Chr5':
                chr5 += 1
        lijst_grootte_chr = [chr1, chr2, chr3, chr4, chr5]  # lijst met de waarde voor elke staaf van de staafdiagram
        hoogte = lijst_grootte_chr
        plt.xlabel('chromosoom nummer')
        plt.ylabel('Hoeveelheid genen')
        plt.title('aantal genen per chromosoom')
        plt.bar(range(len(hoogte)), hoogte, tick_label=chromosomen_nummers)
        plt.show()
    except (ValueError, NameError):
        print('Getal klopt niet kan bijoorbeeld een string/woord zijn ipv een getal/int of naam van var is fout')


class Gui:
    """
    deze functie maakt een gui pop up scherm met alle genen met een drop down menu
    :return: GUI
    """

    def __init__(self):
        self.window = tkinter.Tk()
        self.boven = tkinter.Frame(self.window)
        self.onder = tkinter.Frame(self.window)
        self.rechts = tkinter.Frame(self.window)
        self.id = ['AT1G01010', 'AT1G01020', 'AT1G01030']  # geeft de gen namen voor het dropdown menu
        self.variable = StringVar(self.boven)
        self.variable.set(self.id[0])
        self.dropdown = tkinter.OptionMenu(self.boven, self.variable, *self.id)
        self.dropdown.pack(side='top')
        self.knop2 = tkinter.Button(self.rechts, text='Quit', command=self.window.destroy)  # maakt een stop knop
        self.waarde = tkinter.StringVar()
        self.waarde_label = tkinter.Label(self.boven)
        self.knop2.pack(side='right')
        self.boven.pack()
        self.onder.pack()
        self.rechts.pack()
        self.waarde_label.pack()
        self.labelTest = tkinter.Label(text="", font=('Helvetica', 12), fg='red')
        self.labelTest.pack(side="top")
        tkinter.mainloop()


show = Gui()

main()
