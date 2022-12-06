import fpdf as FPDF



class PDF(FPDF.FPDF):

    def __init__(self,title,author):
        super(PDF,self).__init__()
        self.figcount=1
        FPDF.FPDF.__init__(self,orientation='P',unit='mm',format='A4')
        self.set_font('Arial','',15)
        self.fontsize=15
        self.title=title
        self.add_page()
        self.alias_nb_pages()
        self.set_title(title)
        self.set_author('Remy Zyngfogel')
    def header(self):
        # Arial bold 15
        self.set_font('Arial', 'B', 15)
        # Calculate width of title and position
        w = self.get_string_width(self.title) + 6
        self.set_x((210 - w) / 2)
        # Colors of frame, background and text
        self.set_draw_color(0, 0, 0)
        self.set_fill_color(256, 256, 256)
        self.set_text_color(0, 0, 0)
        # Thickness of frame (1 mm)
        self.set_line_width(.5)
        # Title
        self.cell(w, 9, self.title, 1, 1, 'C', 1)
        # Line break
        self.ln(10)

    def add_image(self,title,filename,legend):
        self.ln()
        self.set_font('Arial', '', 15)
        self.cell(200, 10, txt = title, ln = 1, align = 'c')
        self.image(filename,h=210/2,w=297*(2/3),type="png")
        self.set_font('Arial', '', 9)
        if legend:
            self.cell(200, 5, txt = legend, ln = 1, align = 'c')

    # Page footer
    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Arial italic 8
        self.set_font('Arial', 'I', 8)
        # Page number
        self.cell(0, 10, 'Page ' + str(self.page_no()) + '/{nb}', 0, 0, 'C')


    def add_dict(self,title,message):
        self.ln()
        self.set_font('Arial', '', 15)
        self.cell(200, 10, txt = title, ln = 1, align = 'l')
        self.set_font('Arial', '', 12)

        #self.cell(0,11/self.k, 'Printing line number ' , 0, 1)
        # insert the texts in pdf 
        for i,key in enumerate(message.keys()): 
            sentence='    %s: %s' % (key,str(message[key]))
            self.cell(200, 5, txt = sentence, ln = i+2, align = 'l') 


