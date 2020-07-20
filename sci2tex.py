def getBeforeDigit(n_):
    n = n_.lower()
    i = n.find('.')
    if i == -1:
        i = n.find('e')
        if i == -1:
            return n
    return n[:i]

def getAfterDigitBeforeE(n_):
    n = n_.lower()
    i1 = n.find('.')
    i2 = n.find('e')
    if i1 == -1 and i2 == -1:
        return ''
    if i1 == -1 and i2 != -1:
        return ''
    if i1 != -1 and i2 == -1:
        return n[i1:]
    return n[i1:i2]

def getAfterE(n_):
    n = n_.lower()
    i = n.find('e')
    if i != -1:
        return n[i+1:]
    return ''

def sci2tex(n):
    bfd = getBeforeDigit(n)
    afd = getAfterDigitBeforeE(n)
    afe = getAfterE(n)
    bfd = bfd.lstrip('+').lstrip('0')
    afd = afd.rstrip('0')
    afe = afe.lstrip('+').lstrip('0')
    if afe == '':
      if len(afd) <= 1:
        if len(bfd) == 0:
            return '0'
        return bfd
      return bfd + afd
    if afe.startswith('-'):
        afe = '-' + afe[1:].lstrip('0')
    if len(afd) <= 1:
        return bfd + r'{\times} 10^{' + afe + '}'
    return bfd + afd + r'{\times} 10^{' + afe + '}'
