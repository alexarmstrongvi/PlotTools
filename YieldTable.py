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

class UncFloat :
    precision = 2

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
        value = self.value + other.value
        uncertainty = sqrt(self.uncertainty**2 + other.uncertainty**2)
        return UncFloat(value, uncertainty)
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)
    def __sub__(self, other):
        value = self.value - other.value
        uncertainty = sqrt(self.uncertainty**2 + other.uncertainty**2)
        return UncFloat(value, uncertainty)
    def __mul__(self, other):
        value = self.value * other.value
        rel_unc1 = self.uncertainty /abs(self.value)
        rel_unc2 = other.uncertainty / abs(other.value)
        rel_unc = sqrt(rel_unc1**2 + rel_unc2**2)
        uncertainty = rel_unc * abs(value)
        return UncFloat(value, uncertainty)
    def __div__(self, other):
        value = self.value / other.value
        rel_unc1 = self.uncertainty /abs(self.value)
        rel_unc2 = other.uncertainty / abs(other.value)
        rel_unc = sqrt(rel_unc1**2 + rel_unc2**2)
        uncertainty = rel_unc * abs(value)
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


