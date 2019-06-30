class Label :
    def __init__(self, code, txt=None, root=None, tex=None):
        if not txt:
            txt = code
        if not root:
            root = txt
        if not tex:
            tex = txt

        self.code = code
        self.txt = txt
        self.root = root
        self.tex = py_str_to_tex_str(tex)

    def py_str_to_tex_str(s):
        s = s.replace("\\t","\t").replace("\t","\\t")
        s = s.replace("\\b","\b").replace("\b","\\b")
