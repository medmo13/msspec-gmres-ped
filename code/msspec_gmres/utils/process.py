# vim: set et ts=4 sw=4 fdm=indent:
# coding: utf-8

import re
import os
import textwrap


class Patterns(object):
    col  = '(?P<comment> |C|\*)'
    col += '(?P<label>(?: |\d){1,5})'
    col += '(?P<cont> |\d|&)'


    typ  = '(?P<type>'
    typ += 'BYTE|'
    typ += 'CHARACTER(?:\*\(\*\)|\*\d+)?|'
    typ += 'COMPLEX(?:\*8|\*16|\*32)|'
    typ += 'DOUBLE\s+(?:COMPLEX|PRECISION)|'
    typ += 'INTEGER(?:\*2|\*4|\*8)?|'
    typ += 'LOGICAL(?:\*1|\*2|\*4|\*8)?|'
    typ += 'REAL(?:\*4|\*8|\*16)?|'
    typ += 'AUTOMATIC|STATIC)'

    nam  = '[a-z][a-z0-9_]*'
    dim  = '(?:,?.*?(?::.*?)?)+'
    axs  = '[^:]+(:[^:]+)?'




class BaseInfo(object):
    def __init__(self, **kwargs):
        self._attrs = kwargs.keys()
        for kw,val in kwargs.items():
            setattr(self, '_' + kw, val)
    @property
    def info(self):
        s = '=== {}:\n'.format(self.__class__.__name__)    
        for attr in self._attrs:
            s += '    {}: {}\n'.format(attr, repr(getattr(self, attr)))
        return s

    def __repr__(self):
        return '<{}>'.format(self.__class__.__name__)


class DimensionInfo(BaseInfo):
    def __init__(self, **kwargs):
        opts = {'rank': None, 
                'extents': None, 
                'variable': None}
        opts.update(**kwargs)
        BaseInfo.__init__(self, **opts)

    @property
    def rank(self):
        return len(self.extents)
    @property
    def extents(self):
        return self._extents
    @property
    def variable(self):
        return self._variable
    @variable.setter
    def variable(self, value):
        assert isinstance(value, VariableInfo)
        self._variable = value

    def __str__(self):
        s = ''
        for d in self.extents:
            s += d[0]
            if d[1] is not None:
                s += ':' + d[1]
            s += ','
        s = s.strip(',')
        return s

class VariableInfo(BaseInfo):
    def __init__(self, **kwargs):
        opts = {'name': None, 'type': None, 
                'dimension': None, 'subprogram': None}
        opts.update(**kwargs)
        BaseInfo.__init__(self, **opts)
    @property
    def name(self):
        return self._name
    @property
    def type(self):
        return self._type
    @type.setter
    def type(self, value):
        self._type = value
    @property
    def dimension(self):
        return self._dimension
    @dimension.setter
    def dimension(self, value):
        self._dimension = value
    @property
    def subprogram(self):
        return self._subprogram
    @subprogram.setter
    def subprogram(self, value):
        self._subprogram = value

    def __str__(self):
        s = self.name
        if self.dimension is not None:
            s += '(' + str(self.dimension) + ')'
        return s

class FileInfo(BaseInfo):
    def __init__(self, **kwargs):
        opts = {'filename': None, 'content': None, 'subprograms': None}
        opts.update(**kwargs)
        BaseInfo.__init__(self, **opts)
        self._subprograms = find_subprograms('\n'.join(self.content))
        for sp in self._subprograms:
            sp.file = self
        #self._subprograms = None
    @property
    def filename(self):
        return os.path.abspath(self._filename)
    @property
    def content(self):
        with open(self.filename, 'r') as fd:
            lines = fd.readlines()
            pat = re.compile('     (?:\d|&)\s*(.*)$')
            c = []
            for line in lines:
                line = line.strip('\n')
                m = pat.match(line)
                if m:
                    c[-1] += m.group(1)
                else:
                    c.append(line)
        return c 
    @property
    def subprograms(self):
        return self._subprograms
    def __str__(self):
        return content2str(self.content)

class SubProgramInfo(BaseInfo):
    def __init__(self, **kwargs):
        opts = {'name': None, 'content': None, 'type': None, 'file': None, 
                'l0': None, 'l1': None, 'commons': None}
        opts.update(**kwargs)
        BaseInfo.__init__(self, **opts)

    @property
    def name(self):
        return self._name
    @property
    def content(self):
        return self._file.content[self.l0:self.l1]
    @property
    def type(self):
        return self._type
    @property
    def file(self):
        return self._file
    @property
    def l0(self):
        return self._l0
    @property
    def l1(self):
        return self._l1

    @file.setter
    def file(self, value):
        self._file = value
    
    @property
    def commons(self):
        self._commons = find_commons('\n'.join(self.content))
        for c in self._commons:
            c.subprogram = self
        return self._commons

    def __str__(self):
        return content2str(self.content)

    def __repr__(self):
        s = "{}<{}>".format(self.name, self.type)
        return s


class CommonBlockInfo(BaseInfo):
    def __init__(self, **kwargs):
        opts = {'name': None, 'content': None, 'subprogram': None, 'variables': None}
        opts.update(**kwargs)
        BaseInfo.__init__(self, **opts)

    def __str__(self):
        return content2str(self.content)

    def __repr__(self):
        s = "{}<{:d} variables>".format(self.name, len(self.variables))
        return s
    
    @property
    def subprogram(self):
        return self._subprogram
    @subprogram.setter
    def subprogram(self, value):
        self._subprogram = value
        

    @property
    def name(self):
        return self._name
                
    @property
    def content(self):
        return self._content
    @property
    def variables(self):
        # string to analyse
        s = self.content[0]
        m = re.match("^.*COMMON\s*/{}/\s*(.*)$".format(self.name),s, re.I)
        self._variables = find_variables(m.groups()[0])
        for v in self._variables:
            v.subprogram = self.subprogram

        # If the dimension of a variable is None, try to search if a DIMENSION
        # statement exists in the subprogram
        dim_defs = [] # list of variables defined in a DIMENSION statement
        for line in self.subprogram.content:
            m = re.match("^[^C]\s+DIMENSION\s+(.*)$", line, re.I)
            if m is not None:
                s = m.groups()[0]
                var_list = find_variables(s)
                for v in var_list:
                    dim_defs.append(v)
        dim_defs = Variables(dim_defs)
        # Now for each variable of the common, if there is no dimension, try to find it in
        # the dim_defs list
        for v in self._variables:
            if v.dimension is None:
                dim = dim_defs[v.name]
                if dim is not None:
                    #print(self.subprogram.name)
                    #print(v.name, dim)
                    v.dimension = dim.dimension
                    #exit()

        # if the type of the variable is None, try to find it in a declaration statement
        type_defs = [] # list of variables defined with their type
        for line in self.subprogram.content:
            var_list = find_type(line)
            if var_list is not None:
                for v in var_list:
                    type_defs.append(v)
        type_defs = Variables(type_defs)
        # Now for each variable with no type, try to find it in the type_defs list
        for v in self._variables:
            if v.type is None:
                typ = type_defs[v.name]
                if typ is not None:
                    #print(self.subprogram.name)
                    #print(v.name, typ)
                    v.type = typ.type
                    #exit()
                else:
                    if re.match("^[A-HO-Z].*", v.name, re.I):
                        v.type = "REAL"
                    else:
                        v.type = "INTEGER"
                

        return self._variables

            
class _CommonBlockInfo(BaseInfo):
    def __init__(self, **kwargs):
        self.name = kwargs.get('name', None)
        self.subprogram = kwargs.get('subprogram', None)

    def __init__(self, **kwargs):
        opts = {'name': None, 'content': None, 'type': None, 'file': None, 
                'l0': None, 'l1': None}
        opts.update(**kwargs)
        BaseInfo.__init__(self, **opts)

    def __str__(self):
        s  = "Common block name: {}\n".format(self.name)
        for var in self.variables:
            s += str(var)
        return s
    def __repr__(self):
        s = '{}'.format(self.name)
        return s

    def _find_variables(self):
        content = self.content[0].rstrip()
        
        pat = re.compile('^\s*COMMON\s+/{}/\s+(.*)$'.format(self.name), re.IGNORECASE)
        m = pat.match(content)
        var_loc = m.group(1)

        p0 = re.compile(r'\(.*[^\)]$')
        p1 = re.compile('^[a-zA-Z0-9_\*]*\)$')
        var_list = []
        var_loc_list = var_loc.split(',')
        for i, _ in enumerate(var_loc_list):
            _ = _.strip()
            if i > 0:
                if p0.search(var_list[-1]):
                    var_list[-1] += ',' + _
                else:
                    var_list.append(_)
            else:
                var_list.append(_)

        variables = []
        for var in var_list:
            m = re.match("([a-zA-Z0-9_]+)(\((.*)\))?", var)
            var_name = m.group(1)
            v = VariableInfo(name=var_name, subprogram=self.subprogram)
            
            #dim_loc = m.groups()[-1]
            #if dim_loc is not None:
            #    dim_list = dim_loc.split(',')
            #    v.dimension = dim_list
            variables.append(v)
        return variables
                
    @property
    def content(self):
        pat = re.compile('\s+.*COMMON\s+/{}/.*$'.format(self.name), re.IGNORECASE)
        for line in self.subprogram.content:
            if pat.match(line):
                return [line.rstrip(),]
    @property
    def variables(self):
        variables_list = self._find_variables()
        return Variables(variables_list)




class InfoList(object):
    def __init__(self, elements):
        self._elements = elements
        for element in self._elements:
            setattr(self, element.name, element)
    def __getitem__(self, item):
        if isinstance(item, str):
            for element in self._elements:
                if element.name == item:
                    return element
            return None
        elif isinstance(item, int):
            return self._elements[item]
        else:
            raise NameError('Unable to retrieve item!')
    def __str__(self):
        s  = "nb of {}: {:d}\n".format(self.__class__.__name__, len(self._elements))
        #s += str([_.name for _ in self._elements])
        for _ in self._elements:
            s += _.info
        return s
    def __len__(self):
        return len(self._elements)

class Subprograms(InfoList):
    def __init__(self, *args, **kwargs):
        InfoList.__init__(self, *args, **kwargs)

class Commons(InfoList):
    def __init__(self, *args, **kwargs):
        InfoList.__init__(self, *args, **kwargs)

class Variables(InfoList):
    def __init__(self, *args, **kwargs):
        InfoList.__init__(self, *args, **kwargs)







"""
class Subprograms(BaseInfo):
    def __init__(self, subprograms_list):
        self._subprograms = subprograms_list
        for sp in self._subprograms:
            setattr(self, sp.name, sp)
    def __getitem__(self, item):
        return self._subprograms[item]
    def __str__(self):
        s = "nb of subprograms: {:d}\n".format(len(self._subprograms))
        return s

class Commons(BaseInfo):
    def __init__(self, commons_list):
        self._commons = commons_list
        for cmn in self._commons:
            setattr(self, cmn.name, cmn)
    def __getitem__(self, item):
        return self._commons[item]
    def __str__(self):
        s = "nb of commons: {:d}\n".format(len(self._commons))
        return s

class Variables(BaseInfo):
"""

    


class VariableInfo2(BaseInfo):
    def __init__(self, **kwargs):
        self.name = kwargs.get('name', None)
        #self.type = kwargs.get('type', None)
        self.dimension = kwargs.get('dimension', None)
        self.i = []
        self.j = []
        self.subprogram = kwargs.get('subprogram', None)
        #self._find_type()
        #self._find_dimension()
    def __str__(self):
        s  = "Variable name: {}\n".format(self.name)
        s += "         type: {}\n".format(self.type)
        s += "    dimension: {}\n".format(self.dimension)
        return s

    def _find_implicit(self):
        content = self.subprogram.content
        pat = re.compile('^\s*IMPLICIT\s+', re.IGNORECASE)

    def _find_type(self):
        content = self.subprogram.content
        pat = re.compile('^\s+((?:INTEGER|REAL|DOUBLE PRECISION|COMPLEX|LOGICAL|CHARACTER)\S*).*{}[\(,]?.*$'.format(self.name), re.IGNORECASE)
        for line in content:
            m = pat.match(line)
            print(line)
            if m:
                return m.group(1).strip()
        return None

    def _find_dimension(self):
        content = self.subprogram.content
        dimension = None
        #pat = re.compile('^\s+DIMENSION.*{}\(([^\(])\).*$'.format(self.name),re.IGNORECASE)
        pat = re.compile('^\s+DIMENSION.*{}\((.*?)\).*$'.format(self.name),re.IGNORECASE)
        for line in content:
            m = pat.match(line)
            if m:
                print(line)
                #dimension = m.group(1).strip()
                dimension = m.group(1)
                print(dimension)
        if dimension is not None:
            print('Variable: {}, dimension: {}'.format(self.name, dimension))
    @property
    def type(self):
        t = self._find_type()
        return t

###############################################################################
# UTILITY FUNCTIONS
###############################################################################
    
def splitline_(line, width=72):
    result = []
    i = 0
    j = len(line)
    ll = line[i:j]

    if len(ll) > width:
        s = ''
        for dec in range(8):
            s += '{:d}         '.format(dec)
        print(s)
        print('0123456789'*7+'012')
    #print(line)

    while len(ll) > width:
        breaks = [_.end() for _ in re.finditer('[ ,]', ll)]
        print(ll)
        print(breaks)

        for ij, j in enumerate(breaks):
            tmp = ll[i:j]
            #print(breaks,i,j,ij, tmp, len(tmp))
            if len(tmp) <= width and ij < len(breaks)-1:
                continue
            else:
                _ = ll[i:breaks[ij-1]]
                result.append(_)
                i = len(_)
                if i <= 6:
                    print(j, _)
                    raise NameError('Impossible to cut line at breaks')
                ll = '     &' + ll[i:]
                i = 0
                break
    result.append(ll)
    print(ll)
    return result

def splitline(line, width=72):
    if len(line) == 0:
        return []
    if line[0].upper() == 'C':
        return [line,]
    head = line[:6]
    L = line[6:]     # the working line
    # find the indentation
    m = re.search('^\s*', L)
    indent = L[m.start():m.end()]
    # and define the true width to work with
    W = width - 6 - len(indent)
    
    def rule():
        s = ''
        for dec in range(20):
            s += '{:d}         '.format(dec)
        print(s)
        print('0123456789'*20)

    # find all places to break the line
    breaks = [_.end() for _ in re.finditer('[ ,()+-/*=:]', L)]

    # split at breaks
    indices = [0,] + breaks + [len(L),]
    indices = zip(indices[:-1], indices[1:])
    splitted_line = [L[a:b] for a,b in indices]

    # iterate over each element and add it to the previous one if length is < max
    chain = [splitted_line[0],]
    for element in splitted_line[1:]:
        l1 = len(chain[-1])
        l2 = len(element)
        if l1+l2 < W:
            chain[-1] = chain[-1] + element
        else:
            chain.append(element)
    # restore the head of the line
    chain[0] = head + chain[0]
    # add the & symbol
    for i in range(1,len(chain)):
        chain[i] = "     &" + indent + chain[i]
    
    # final check
    for element in chain:
        if len(element) > width:
            rule()
            print(f"{line}")
            print(f"breaks at = {breaks}")
            print(chain)
            rule()
            print(element)
            raise NameError(f"Unable to split line!")
    
        

    return chain



def content2str(content):
    new_content = []
    for index, line in enumerate(content):
        #print(f'{index:>5d}#{line}')
        multilines = splitline(line.rstrip(), width=72)
        #print(multilines)
        new_content += multilines
        
    return '\n'.join(new_content)   


def split_at_comma(string):
    """
    """
    # remove all spaces from the string
    line0 = string.replace(' ', '')
    line  = line0
    # define some patterns
    pat0 = re.compile('\(([^\(\)]*)\)', re.I)

    # remove nested blocks in ()'s and replace them with
    # hash signs (#) to make the treatment easier
    while True:
        M = list(pat0.finditer(line))
        if len(M) == 0: break
        for m in M:
            i,j = m.start(), m.end()
            line = line[:i] + '#'*(j-i) + line[j:]

    # now get indices of ','
    indices = [_.start() for _ in re.finditer(',', line)]
    indices = zip([-1,] + indices, indices + [len(line),])
    # create the list of
    elements = [line0[i+1:j] for i,j in indices]
    return elements

def find_dimension(string):
    """ 
    Finds the components of a dimension declaration.

    :param string: The argument of a dimension declaration
    :type string: str
    :return: A DimensionInfo object.
    :rtype: DimensionInfo

    :Example:

    >>> dim = find_dimension('I,J,-3:2')
    >>> print(dim)
    >>> (3, ())


    """
    # define some patterns
    pat0 = re.compile('([^:]+):?([^:]+)?', re.I)

    # create the list of axes
    axl = split_at_comma(string)
    # get the extents
    extents = []
    for ax in axl:
        m = pat0.match(ax)
        extents.append(m.groups())

    return DimensionInfo(extents=extents)

def find_variables(string):
    """ 
    Finds the name and dimension of variables in a comma separated
    list of variables.

    :param string: The comma separated variables declaration
    :type string: str
    :return: A Variables object.
    :rtype: Variables

    :Example:

    >>> variables = find_variables('ONE, TWO(3,3)')
    >>> print(variables)
    >>> nb of Variables: 2
    >>> === VariableInfo:
    >>>     name: 'ONE'
    >>>     type: None
    >>>     dimension: None
    >>>     subprogram: None
    >>> === VariableInfo:
    >>>     name: 'TWO'
    >>>     type: None
    >>>     dimension: <DimensionInfo>
    >>>     subprogram: None

    """
     # create the list of variables
    var_list = []
    variables = split_at_comma(string)
    pat0 = re.compile('({})(?:\((.*)\))?'.format(Patterns.nam), re.I)
    for var in variables:
        # extract the variable's name and dimension if any
        m = pat0.match(var)
        name = m.group(1)
        if m.group(2) is not None:
            dimension = find_dimension(m.group(2))
        else:
            dimension = None

        variable = VariableInfo(name=name, dimension=dimension)
        if isinstance(dimension, DimensionInfo):
            dimension.variable = variable
        var_list.append(variable)

    return Variables(var_list)

def find_subprograms(string):
    lines = string.split('\n')
    subprograms = []
    for iline, line in enumerate(lines):
        patterns = [('SUBROUTINE', re.compile("\s*SUBROUTINE\s+(\w+)\(?.*")),
                    ('FUNCTION', re.compile("\s*.*FUNCTION\s+(\w+)\(?.*")),
                    ('PROGRAM', re.compile("\s*PROGRAM\s+(\w+).*"))]
        for t, pat in patterns:
            m = pat.match(line)
            if m is not None:
                subprog = SubProgramInfo(type=t, 
                                         name=m.group(1),
                                         l0=iline)
                subprograms.append(subprog)

    for i, subprog in enumerate(subprograms):
        if i < len(subprograms) - 1:
            subprog._l1 = subprograms[i+1].l0
        else:
            subprog._l1 = -1

    return Subprograms(subprograms)

def find_commons(string):
    pat = re.compile("^\s+COMMON\s*/([a-zA-Z0-9_]+)/(.*)$")
    commons = []
    for line in string.split('\n'):
        # extract the name of the common block
        m = pat.match(line)
        if m is not None:
            # name
            name = m.group(1)
            c = CommonBlockInfo(name=name, content=[line,])
            commons.append(c)
    return Commons(commons)


def find_names(string):
    """
    Find the names in expression ie remove (,),+,-,*,/,**,=
    """
    m = re.findall('[A-Z][A-Z0-9_]*', string, re.I)
    return set(m)


def find_type(string):
    """
    return a Variables object if string is a type declaration
    """
    # get out if string is a comment
    if string.upper().startswith('C'):
        return None

    pat  = "^\s+"
    pat += "(BYTE|"
    pat += "CHARACTER(?:\*[0-9]+)?|CHARACTER\*\(\*\)|"
    pat += "COMPLEX(?:\*(?:8|16|32))?|"
    pat += "DOUBLE COMPLEX|"
    pat += "DOUBLE PRECISION|"
    pat += "INTEGER(?:\*(?:2|4|8))?|"
    pat += "LOGICAL(?:\*(?:1|2|4|8))?|"
    pat += "REAL(?:\*(?:4|8|16))?|"
    pat += "AUTOMATIC|"
    pat += "STATIC)"
    pat += "\s+(.*)$"
    m = re.match(pat, string, re.I)
    if m is not None:
        if re.search("IMPLICIT", string):
            return None
        else:
            var_list = find_variables(m.groups()[1])
            for var in var_list:
                var.type = m.groups()[0]
            return var_list

def find_dim(string):
    """
    return a Variables object if string is a dimension declaration
    """
    pat  = "^\s+DIMENSION\s+(.*)$"
    m = re.match(pat, string, re.I)
    if m is not None:
        var_list = find_variables(m.groups()[0])
        return var_list


def write_modules(infile):
    fi = FileInfo(filename=infile)

    # Get all the common blocks defined in the source file
    all_commons = []
    for sp in fi.subprograms:
        for c in sp.commons:
            if c.name not in [_.name for _ in all_commons]:
                all_commons.append(c)

    all_commons = Commons(all_commons)

    # a function to create a module fortran code from a CommonBlockInfo object
    def common2module(cbi):
        variables = cbi.variables
        alloc_args = set()
        module_name = cbi.name.upper() + "_MOD"       


        for variable in variables:
            if variable.dimension is not None:
                dim_vars = find_names(str(variable.dimension))
                alloc_args.update(dim_vars)

 
        s  = f"MODULE {module_name}\n"
        s +=  "  IMPLICIT NONE\n"

        # for each variable whose type is defined explicitely
        for variable in variables:
            s += f"  {variable.type}"
            dimension = variable.dimension
            if dimension is not None:
                if dimension.rank > 0:
                    s += f", ALLOCATABLE, DIMENSION(:" + ",:" * (variable.dimension.rank-1) + ")"
            s += f" :: {variable.name}\n"

        s +=  "CONTAINS\n"
        #s += f"  SUBROUTINE ALLOC_{cbi.name.upper()}({','.join(alloc_args)})\n"
        #s +=  "    IMPLICIT INTEGER (A-Z)\n"
        s += f"  SUBROUTINE ALLOC_{cbi.name.upper()}()\n"
        s += f"    USE DIM_MOD\n"
        # for each variable with a defined dimension
        for variable in variables:
            if variable.dimension is not None:
                s += f"    IF (ALLOCATED({variable.name})) THEN\n"
                s += f"      DEALLOCATE({variable.name})\n"
                s += f"    ENDIF\n"
                s += f"    ALLOCATE({variable})\n"
        
        s += f"  END SUBROUTINE ALLOC_{cbi.name.upper()}\n"
        s += f"END MODULE {module_name}\n"

        # indentation
        s = textwrap.indent(s, prefix="      ")
        # split in too long
        content = s.split('\n')
        s = content2str(content)
        return s

    # write the modules.f file
    with open("modules.f", "w") as fd:
        for c in all_commons:
            fd.write("C" + "="*71 + "\n")
            s = common2module(c)
            fd.write(s)
            fd.write("\n"*2)




if __name__ == "__main__":
    infile = 'spec.f.orig'

    # write the modules.f file
    write_modules(infile)
    exit()

    fi = FileInfo(filename=infile)

    # Get all the common blocks defined in the source file
    all_commons = []
    for sp in fi.subprograms:
        for c in sp.commons:
            if c.name not in [_.name for _ in all_commons]:
                all_commons.append(c)

    all_commons = Commons(all_commons)
    content = fi.content

    # write the allocation.f file
    dim_vars = [
        "NATP_M", 
        "NATCLU_M_", 
        "NAT_EQ_M", 
        "N_CL_L_M",
        "NE_M", 
        "NL_M",
        "LI_M", 
        "NEMET_M",
        "NO_ST_M",
        "NDIF_M",
        "NSO_M",
        "NTEMP_M",
        "NODES_EX_M",
        "NSPIN_M",
        "NTH_M", 
        "NPH_M",
        "NDIM_M",
        "N_TILT_M",
        "N_ORD_M",
        "NPATH_M",
        "NGR_M"]

    s  = f"SUBROUTINE ALLOCATION({', '.join([v+'_' for v in dim_vars])})\n"
    s += f"  USE DIM_MOD\n"
    s += f"  IMPLICIT INTEGER (A-Z)\n"
    for v in dim_vars:
        s += f"  {v} = {v+'_'}\n"
    s += f"  CALL INIT_DIM()\n"
    s += f"END SUBROUTINE ALLOCATION\n"

    # indentation
    s = textwrap.indent(s, prefix="      ")
    # split in too long
    s = content2str(s.split('\n'))
    
    with open("allocation.f", "w") as fd:
        fd.write(s)



    #exit()


    # remove type definitions for variables that are in commons
    for sp in fi.subprograms:
        # get the list all all variables in all commons in this subprogram
        vlist = []
        for c in sp.commons:
            for v in c.variables:
                vlist.append(v.name)
        for iline, line in enumerate(sp.content):
            print(f"{sp.l0+iline:05d}: {line}")
            newline = ''
            # comment INCLUDE statement
            if re.search('INCLUDE.*spec.inc', line, re.I):
                newline = "C" + line
            # replace the commons by USE statements
            m = re.match("^.*COMMON\s+/(.*)/.*$", line, re.I)
            if m is not None:
                cmn_name = m.groups()[0]
                # Here we test if the common variables are the same than the module
                cmn_variables = sp.commons[cmn_name].variables
                s = f"      USE {cmn_name.upper()}_MOD "
                modifications = []
                for i, v in enumerate(cmn_variables):
                    original = all_commons[cmn_name].variables[i].name
                    if v.name != original:
                        modifications.append(f"{v.name} => {original}")
                s += ", ".join(modifications)
                newline = s

            # Remove type declaration for variables that are now in modules
            allv = find_type(line)
            newallv = []
            line_ = line
            if allv is not None:
                # Here the line is a declaration statement
                # remove every variables that are also in vlist
                for v in allv:
                    if v.name not in vlist:
                        # keep this variable in the list
                        newallv.append(str(v))
                        
                # if there is no change
                if len(allv) == len(newallv):
                    line_ = ""
                
                # if no more variables are defined, remove the line
                elif len(newallv) == 0:
                    line_ = "C"
                else:
                    line_ = "      " + v.type + " " + ",".join(newallv)
                
                newline = line_

            # Remove dimension declaration for variables that are now in modules
            allv = find_dim(line)
            newallv = []
            line_ = line
            if allv is not None:
                # Here the line is a dimension statement
                # remove every variables that are also in vlist
                for v in allv:
                    if v.name not in vlist:
                        # keep this variable in the list
                        #if v.dimension is not None:
                        newallv.append(str(v))
                        
                # if there is no change
                if len(allv) == len(newallv):
                    line_ = ""
                # if no more variables are defined, remove the line
                elif len(newallv) == 0:
                    line_ = "C"
                else:
                    line_ = "      DIMENSION " + ",".join(newallv)
                
                newline = line_

            if newline != '':
                print(sp.l0, iline, sp.l0+iline)
                content[sp.l0+iline] = newline
                print(f">>> {newline}")
                    



    # rewrite the file
    with open("new.f", "w") as fd:
        fd.write(content2str(content))
     





