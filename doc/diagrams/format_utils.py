def uml_formatter(class_name,
                  methods=[],
                  members=[]):
    label = '{%s|\l' % class_name
    for method in methods:
        label += ' %s\l\\n' % method
    label += '|\l'
    for member in members:
        label += ' %s\l\\n' % member

    label += '}'
    label = label.replace('<', '\<')
    label = label.replace('>', '\>')
    label = label.replace(':', '\:')
    return label