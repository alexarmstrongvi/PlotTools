"""
================================================================================
Class for storing and manipulating yields

Author:
    Alex Armstrong <alarmstr@cern.ch>

License:
    Copyright: (C) <May 20th, 2018>; University of California, Irvine
================================================================================
"""

# General python
import sys, os
import re
from numbers import Number
from math import sqrt
from collections import OrderedDict, defaultdict
from pdb import set_trace
from tabulate import tabulate

class UncFloat :
    precision = 2
    #Assume all errors are uncorrelated

    def __init__(self, value = 0 , uncertainty = 0):
        if isinstance(value, Number):
            pass # Expected
        elif isinstance(value, basestring):
            value, uncertainty = parse_uncertainty_string(value)
        else:
            print "WARNING : Unexpected input type"

        self.value = float(value)
        self.uncertainty = float(uncertainty)

    def __add__(self, other):
        if isinstance(other, Number):
            value = self.value + other 
            uncertainty = self.uncertainty
        else:
            value = self.value + other.value
            uncertainty = sqrt(self.uncertainty**2 + other.uncertainty**2)
        return UncFloat(value, uncertainty)
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)
    def __sub__(self, other):
        if isinstance(other, Number):
            value = self.value - other 
            uncertainty = self.uncertainty
        else:
            value = self.value - other.value
            uncertainty = sqrt(self.uncertainty**2 + other.uncertainty**2)
        return UncFloat(value, uncertainty)
    def __mul__(self, other):
        if isinstance(other, Number):
            value = self.value * other 
            uncertainty = self.uncertainty * other
        else:
            # Add relative error in quadrature
            # Solve for relative uncertainty to avoid problem from nominal = 0
            # This would lead to undefined relative error (e.g. uncX/X = uncX/0 = undefined)
            uncertainty = sqrt((self.value * other.uncertainty)**2 + (other.value * self.uncertainty)**2)
        return UncFloat(value, uncertainty)
    def __rdiv__(self, other):
        if isinstance(other, Number):
            # Assume errors are extremes
            # Use largest uncertainty up/down
            value = other / self.value
            up = self.value + self.uncertainty
            dn = self.value - self.uncertainty
            if dn < 0 and 0 < up:
                print "WARNING :: Dividing by number consistent with zero:", self
                print "INFO :: Setting error to infinity"
                uncertainty = float("inf")
            else:
                unc_up = abs((value) - (other / up))
                unc_dn = abs((value) - (other / dn))
                uncertainty = max(unc_up, unc_dn)
            return UncFloat(value, uncertainty)
        else:
            raise TypeError, "Undefined division of UncFloat by %s " % str(type(other))
    def __div__(self, other):
        if isinstance(other, Number):
            value = self.value / other 
            uncertainty = self.uncertainty / other
        else:
            value = self.value / other.value
            dn = other.value - other.uncertainty
            up = other.value + other.uncertainty
            if dn < 0 and 0 < up:
                print "WARNING :: Dividing by number consistent with zero:", other
                print "INFO :: Setting error to infinity"
                uncertainty = float("inf")
            else:
                # Add relative error in quadrature
                # Solve for relative uncertainty to avoid problem from nominal = 0 in numerator
                # This would lead to undefined relative error (e.g. uncX/X = uncX/0 = undefined)
                x1 = self.uncertainty * (1.0 / other.value)
                x2 = other.uncertainty * (-self.value / other.value**2)
                uncertainty = sqrt(x1**2 + x2**2)
        return UncFloat(value, uncertainty)
    def __lt__(self, other) :
        return self.value < other.value
    def __le__(self, other) :
        return self.value <= other.value
    def __eq__(self, other) :
        return (isinstance(other, self.__class__)
            and self.value == other.value
            and self.uncertainty == other.uncertainty)
    def __ne__(self, other) :
        return not self.value == other.value
    def __gt__(self, other) :
        return self.value > other.value
    def __ge__(self, other) :
        return self.value >= other.value

    def __str__(self):
        return "%.*f +/- %.*f"%(self.precision, self.value,
                                self.precision, self.uncertainty)

    def parse_uncertainty_string(string):
        if '+/-' not in string:
            print "WARNING : Unrecognized uncertainty format:", string
            return 0, 0

        nums = [n.strip() for n in string.split('+/-')]
        nums = (float(n) for n in nums if self.is_number(n))
        if len(nums) != 2:
            print "WARNING : Unrecognized uncertainty format:", string
            return 0, 0

        return nums[0], nums[1]

    def is_number(string):
        try:
            float(string)
            return True
        except ValueError:
            return False

class TableLabel:
    def __init__(self, name, display="", latex=""):
        self.name = name
        self.display = display
        self.latex = latex

        if not display: self.display = self.name
        if not latex: self.latex = self.display 

class YieldTbl:
    precision = 2
    write_to_latex = False

    def __init__(self):
        
        self.row_names = []
        self.row_displaynames = []
        self.row_latexnames = []
        self.row_mc_flags = []
        self.row_signal_flags = []
        
        self.col_names = []
        self.col_displaynames = []
        self.col_latexnames = []
        self.common_column_prefix = ""
        
        self.table = []

        self.column_formulas = []
        self.row_formulas = []

    def add_row(self, row_name, row_displayname = "", row_latexname="", mc=False, signal=False):
        display = row_displayname if row_displayname else row_name
        latex = row_latexname if row_latexname else display

        self.row_names.append(row_name)
        self.row_displaynames.append(display)
        self.row_latexnames.append(latex)
        self.row_mc_flags.append(mc)
        self.row_signal_flags.append(signal)

        n_cols = len(self.table[0]) if len(self.table) else 1
        row = [UncFloat(0,0) for x in range(n_cols)] 
        self.table.append(row) 
        

    def add_column(self, col_name, col_displayname = "", col_latexname=""):
        display = col_displayname if col_displayname else col_name
        latex = col_latexname if col_latexname else display

        self.col_names.append(col_name)
        self.col_displaynames.append(display)
        self.col_latexnames.append(latex)

        for row in self.table:
            row.append(UncFloat(0,0))        

    def add_entry(self, row_name, col_name, val,
                  error=None,
                  row_displayname="",
                  row_latexname="",
                  col_displayname="",
                  col_latexname="",
                  mc = False,
                  signal = False,
                  ):
        if not row_displayname: row_displayname = row_name
        if not row_latexname: row_latexname = row_displayname
        if not col_displayname: col_displayname = col_name
        if not col_latexname: col_latexname = col_displayname
        
        if row_name not in self.row_names:
            self.add_row(row_name, row_displayname, row_latexname, mc = mc, signal=signal)
            if col_name not in self.col_names:
                self.col_names.append(col_name)
                self.col_displaynames.append(col_displayname)
                self.col_latexnames.append(col_latexname)
        
        if col_name not in self.col_names:
            self.add_column(col_name, col_displayname, col_latexname)

        
        try:
            row_idx = self.row_names.index(row_name)
            col_idx = self.col_names.index(col_name)
            self.table[row_idx][col_idx] = UncFloat(val, error)
        except ValueError:
            print "ERROR: row %s and column %s are not found in this yield table" % (row_name, col_name)
            print "INFO: Current rows:", self.row_names 
            print "INFO: Current columns:", self.col_names 

    def get(row_name, col_name):
        try:
            row_idx = self.row_names.index(row_name)
            col_idx = self.col_names.index(col_name)
            return table[row_idx][col_idx]
        except ValueError:
            print "ERROR: row %s and column %s are not found in this yield table" % (row_name, col_name)
            print "INFO: Current rows:", self.row_names 
            print "INFO: Current columns:", self.col_names 

        
    def get_column(self, col_name):
        col_idx = self.col_names.index(col_name)
        column = []
        for row in self.table:
            column.append(row[col_idx])
        return column

    def add_column_formula(self, name, formula, displayname='', latexname=''):
        label = TableLabel(name, displayname, latexname)
        self.column_formulas.append((label, formula))
        
    def add_row_formula(self, name, formula, displayname='', latexname=''):
        label = TableLabel(name, displayname, latexname)
        self.row_formulas.append((label, formula))
    
    def apply_row_formulas(self):
        new_rows = []
        for formula_pair in self.row_formulas:
            # Apply formula
            tmp_row = self.apply_row_formula(formula_pair[1], formula_pair[0])
            
            # Skip column if all values are NaN
            if all(x != x for x in tmp_row): tmp_row = []
            
            new_rows.append(tmp_row)
        
        # Add formulas to original table
        for r_idx, row in enumerate(new_rows):
            if not row: continue
            for c_idx, entry in enumerate(row):
                if isinstance(entry, float):
                    value, error = entry, 0
                else:
                    value, error = entry.value, entry.uncertainty

                row_names = self.row_formulas[r_idx][0]
                mc = True if row_names.display == 'mc_total' else False
                signal = False
                self.add_entry(row_name=row_names.name,
                               col_name=self.col_names[c_idx],
                               val=value,
                               error=error,
                               row_displayname = row_names.display,
                               row_latexname = row_names.latex,
                               col_displayname = self.col_displaynames[c_idx],
                               col_latexname = self.col_latexnames[c_idx],
                               mc = mc,
                               signal = signal
                               ) 
    
    def apply_row_formula(self, formula, labels):
        # Get rows from formula
        rows = self.extract_elements_from_formula(formula)
        rows = list(set(rows))

        # Replace formula keywords
        for row in rows:
            replace_str = self.replace_formula_keywords(row)
            formula = re.sub(r"\b%s\b"%row, replace_str, formula)
        
        # Replace names in formula with values
        rows = self.extract_elements_from_formula(formula)
        rows = list(set(rows))
        for row in rows:
            try:
                idx = self.row_names.index(row)
                replace_str = "col[%d]" % idx
                formula = re.sub(r"\b%s\b"%row, replace_str, formula)
            except ValueError:
                print "ERROR :: %s not found in row names" % row
                print "INFO :: Row Names:", self.row_names
                formula = "float('NaN')"
        
        new_row = []
        #Idx label must agree with label used in formula string
        for col_name in self.col_names:
            col = self.get_column(col_name)
            if '()' in formula:
                print "ERROR :: Empty term in formula:", formula
                formula = 'float("NaN")'
            try:
                result = eval(formula)
            except ZeroDivisionError:
                result = float('inf')
            except NameError, err:
                print "ERROR :: ", err
                print "INFO :: Row names = ", self.row_names
                result = float('NaN')
            except TypeError, err:
                print "ERROR :: ", err
                print "INFO :: Formula = ", formula 
                result = float('NaN')

            new_row.append(result)
        
        return new_row

    def apply_column_formulas(self):
        new_columns = []
        for formula_pair in self.column_formulas:
            # Apply formula
            tmp_col = self.apply_column_formula(formula_pair[1], formula_pair[0])

            # Skip column if all values are NaN
            if all(x[0] != x[0] for x in tmp_col):
                continue

            # Store new column
            if not new_columns: new_columns = tmp_col
            else:
                for entry, row in zip(tmp_col, new_columns):
                    row.append(entry[0])
        
        # Add formulas to original table
        for r_idx, row in enumerate(new_columns):
            for f_idx, entry in enumerate(row):
                if isinstance(entry, float):
                    value, error = entry, 0
                else:
                    value, error = entry.value, entry.uncertainty

                col_names = self.column_formulas[f_idx][0]
                self.add_entry(row_name=self.row_names[r_idx],
                               col_name=col_names.name,
                               val=value,
                               error=error,
                               row_displayname = self.row_displaynames[r_idx],
                               row_latexname = self.row_latexnames[r_idx],
                               col_displayname = col_names.display,
                               col_latexname = col_names.latex 
                               ) 
    def n_rows(self):
        n1 = len(self.row_names)
        n2 = len(self.row_displaynames)
        n3 = len(self.row_latexnames)
        n4 = len(self.row_mc_flags)
        n5 = len(self.row_signal_flags)
        n6 = len(self.table)
        if n1 == n2 == n3 == n4 == n5:
            return n1
        else:
            print "WARNING :: lost track of row count"
            print "INFO :: ", n1, n2, n3, n4, n5
            return 0

    def mc_sample_names(self, no_signal=True):
        mc_samples = []
        for idx in range(self.n_rows()):
            sig_skip = no_signal and self.row_signal_flags[idx]
            if self.row_mc_flags[idx] and not sig_skip:
                mc_samples.append(self.row_names[idx])
        return mc_samples

    def signal_sample_names(self):
        sig_samples = []
        for row_name, sig_flag in zip(self.row_names, self.row_signal_flags):
            if sig_flag:
                sig_samples.append(row_name)
        return sig_samples
    
    def replace_formula_keywords(self, keyword):
        if not self.common_column_prefix:
            self.common_column_prefix = os.path.commonprefix(self.col_names).strip()
        if keyword == "MC":
            keyword = "(" + "+".join(self.mc_sample_names()) + ")"
        elif keyword == "SIGNAL":
            keyword = "(" + "+".join(self.signal_sample_names()) + ")"
        elif keyword == "BASE_REG":
            keyword = self.common_column_prefix
        elif keyword not in self.col_names:
            new_kw = self.common_column_prefix + "_" + keyword
            if new_kw in self.col_names: 
                keyword = new_kw 
            else:
                # Error raised elsewhere
                pass

        return keyword

    def extract_elements_from_formula(self, formula):
        elements = re.sub("(\+|\-|\)|\(|\*|\.|[0-9]|\/)"," ", formula)
        elements = [c.strip() for c in elements.split()]
        assert all(s.replace("_","").isalpha() for s in elements), (
            "ERROR :: Unacceptable formula format %s -> %s" % (formula, elements))
        return elements

    def apply_column_formula(self, formula, labels):
        # Get columns from formula
        columns = self.extract_elements_from_formula(formula)
        columns = list(set(columns))

        # Replace formula keywords
        for col in columns:
            replace_str = self.replace_formula_keywords(col)
            formula = re.sub(r"\b%s\b"%col, replace_str, formula)
        
        # Replace names in formula with values
        columns = self.extract_elements_from_formula(formula)
        columns = list(set(columns))
        for col in columns:
            try:
                idx = self.col_names.index(col)
                replace_str = "row[%d]" % idx
                formula = re.sub(r"\b%s\b"%col, replace_str, formula)
            except ValueError:
                print "WARNING :: Formula column not stored:", col
                print "INFO :: Col names = ", self.col_names
                formula = "float('NaN')"
        
        new_column = []
        for row in self.table:
            try:
                result = eval(formula)
            except ZeroDivisionError:
                result = float('inf')
            new_column.append([result])
        
        return new_column

    def Print(self, latex=False, mc_data_fmt=False):
        print self.print_str(latex, mc_data_fmt)
    
    def save_table(self, save_path, latex=False, mc_data_fmt=False):
        with open(save_path, "w") as ofile:
            ofile.write(self.print_str(latex, mc_data_fmt))

    def print_str(self, latex=False, mc_data_fmt=False):
        if latex:
            col_names = self.col_latexnames
            row_names = self.row_latexnames
            table_format = 'latex_raw'
            str_f = lambda x : str(x).replace("+/-","$\pm$")
        else:
            col_names = self.col_names
            row_names = self.row_names
            table_format = 'psql'
            str_f = str     
            
        print_table = []
        for idx, row in enumerate(self.table):
            print_table.append([row_names[idx]] + map(str_f, row))

        print_str = tabulate(print_table, headers=col_names, tablefmt=table_format)
        return self.format_print_str(print_str, latex, mc_data_fmt)

    def format_print_str(self, print_str, latex=False, mc_data_fmt=False):
        # Arrange rows
        if mc_data_fmt:
            print_lst = print_str.split("\n")
            if latex:
                end_tabular = print_lst.pop()
                hline = "\hline"
                n_header_rows = 4
            else:
                hline = print_lst[0]
                n_header_rows = 3
            moved_lines = 0
            
            for idx in range(self.n_rows()):
                if self.row_mc_flags[idx] and not self.row_signal_flags[idx]: continue
                line_to_move = idx + n_header_rows - moved_lines
                print_lst.append(print_lst.pop(line_to_move))
                moved_lines += 1

            print_lst.append(hline)
            if latex: 
                print_lst.insert(4,hline)
                print_lst.append(end_tabular)
            print_str = "\n".join(print_lst)
        return print_str




class YieldTable :
    precision = 2

    def __init__(self):
        self.mc = OrderedDict()
        self.data = {}
        self.signals = {}
        self.formulas = {}
        self.region = ""
        self.variable = ""
        self.partitions = [] 
    #TODO: replace data, signals, and mc with a single yields. Turn data, mc, and signals into properties
    # that add the value to yields and flags the type of sample as data, mc, or signal

    def get_mc_yields(self):
        return sum([v for _, v in self.mc.iteritems()])
    
    def all_yields(self):
        all_yields = self.mc.copy()
        all_yields['MC'] = self.get_mc_yields()
        all_yields.update(self.data)
        all_yields.update(self.signals)
        return all_yields

    def all_names(self):
        all_names = []
        all_names += [k for k in self.mc]
        all_names += [k for k in self.signals]
        all_names += [k for k in self.data]
        all_names += [k for k in self.formulas]
        return all_names

    def get_max_mc_yield(self):
        yields = [v for _, v in self.mc.iteritems()]  
        vals = [v.value for v in yields]
        idx = vals.index(max(vals))
        return yields[idx]
         

    def __eq__(self, other):
        ''' Check if all yields, both value and uncertainty, are identical'''
        if not isinstance(other, self.__class__): return false

        y1 = self.all_yields()
        y2 = other.all_yields()
        for (name1, yld1), (name2, yld2) in zip(y1.iteritems(),y2.iteritems()):
            if name1 != name2 or str(yld1) != str(yld2):
                return False
        return True
    
    def __add__(self, other):
        for key in data: self.data[key] += other.data[key]
        for key in mc: self.mc[key] += other.mc[key]
        for key in signals: self.signals[key] += other.signals[key]
            

        value = self.value + other.value
        uncertainty = sqrt(self.uncertainty**2 + other.uncertainty**2)
        return UncFloat(value, uncertainty)
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def is_empty(self):
        return not (self.mc or self.data or self.signals)

    def partition_string(self, name = None, yld = None,
            data = False, mc = False, signals = False, 
            total_mc = False, formulas = False, region = False):
        assert sum([data, mc, signals, total_mc, formulas, region]) == 1, (
                "ERROR :: Must set one and only one boolean")
        assert (yld != None) or (region or formulas)
        assert name or total_mc or region
        partition_str = ""
        remainder = yld
        for ii, part in enumerate(self.partitions):
            if mc and part.mc: 
                val = part.mc[name].value
            elif data and part.data: 
                val = part.data[name].value
            elif signals and part.signals: 
                val = part.signals[name].value
            elif total_mc and part.mc:
                val = part.get_mc_yields().value
            elif region and part.region:
                val = part.region
            elif formulas and part.formulas:
                formula = self.formulas[name]
                val = part.apply_formula(formula)
            else:
                return ""
            op = "+" if ii else "  ="
            if isinstance(val, str):
                partition_str += "%s  %s  "%(op, val)
            else:
                partition_str += "%s  %.*f  "%(op, self.precision, val)
                remainder -= val
        if yld and remainder != 0:
            partition_str += "%s  %.*f  "%(op, self.precision, remainder)
        if self.partitions and region:
            partition_str += "%s  %s  "%(op, 'Other')
        return partition_str
    
    def Print(self):
        mc_total = self.get_mc_yields() 

        formula_values = {}
        for key, formula in self.formulas.iteritems():
            formula_values[key] = self.apply_formula(formula)

        # Get formatting settings
        longest_name = max(self.all_names(), key=len)
        space = len(longest_name) + 2
        longest_mc_yield = self.get_max_mc_yield()
        space2 = len(str(longest_mc_yield)) 

        # Print Table
        print "==============  Yield Table =============="
        print " Variable(s): ", self.variable
        print " Region: ", self.region,
        print self.partition_string(region=True) if self.partitions else ""
        if len(self.mc):
            print "------------------------------------------"
            print "Backgrounds:"
            for name, yield_value in self.mc.iteritems():
                print "%*s : %*s"%(space, name, space2, yield_value),
                print self.partition_string(name=name, yld=yield_value.value, mc=True) if self.partitions else ""
        if len(self.data):
            print "------------------------------------------"
            print "%*s : %*s"%(space, "MC", space2, mc_total),
            print self.partition_string(yld=mc_total.value, total_mc=True) if self.partitions else ""
            for name, yield_value in self.data.iteritems():
                print "%*s : %*s"%(space, name, space2, yield_value.value),
                print self.partition_string(name=name, yld=yield_value.value, data=True) if self.partitions else ""
        if len(self.signals):
            print "------------------------------------------"
            print "Signal:"
            for name, yield_value in self.signals.iteritems():
                print "%*s : %*s"%(space, name, space2, yield_value),
                print self.partition_string(name=name, yld=yield_value.value, signals=True) if self.partitions else ""
        if len(self.formulas):
            print "------------------------------------------"
            for name, value in formula_values.iteritems():
                print "%*s : %*.*f"%(space, name, space2, self.precision, value)
                # TODO: Fix - sampel names are different in different partitions
                #print self.partition_string(name=name, formulas=True) if self.partitions else ""


        print "=========================================="

    def apply_formula(self, formula, no_uncertainty=True):
        # Get samples from formula
        samples = re.sub("(\+|\-|\)|\(|\*|\.|[0-9]|\/)"," ", formula)
        samples = [s.strip() for s in samples.split()]
        assert all(s.replace("_","").isalpha() for s in samples), (
            "ERROR :: Unacceptable formula format %s -> %s" % (formula, samples))


        # Replace names in formula with values
        all_yields = self.all_yields()
        assert len(all_yields) == len(self.mc) + len(self.data) + len(self.signals) + 1, (
            "ERROR (YieldTable) :: Overlapping key values")
        samples = list(set(samples))
        samples.sort(key=len, reverse=True)

        for sample_name in samples:
            if sample_name not in all_yields:
                print "WARNING :: Formula sample not stored:", sample_name
                return float("NaN")
            if no_uncertainty:
                replace_str = "all_yields['%s'].value"%sample_name
            else:
                replace_str = "all_yields['%s']"%sample_name
            formula = re.sub(r"\b%s\b"%sample_name, replace_str, formula)

        # Evaluate formula
        try:
            result = eval(formula)
        except ZeroDivisionError:
            result = float('inf')
        finally:
            return result

    def reset(self):
        self.mc.clear()
        self.data.clear()
        self.signals.clear()
        self.region = ""
        self.variable = ""


