def text2file(text, filename):
    fd = open(filename, 'w')
    fd.write(text)
    fd.close()

def file2text(filename):
    content = []
    with open(filename) as fd:
        for l in fd:
            content.append(l.strip())
    return content
