

import requests
import re
import json


# First kind of page: one page, each trimester has it own section
# This only works if all names are presents.
# Generated files may require manual corrections...
titles = [
    ('fr', 'Calendrier_r%C3%A9publicain', [29, 30, 31, 32], ".*FrRepCalLine.*\\[\\[(.*)\\]\\].*"),
    ('en', 'French_Republican_calendar', [10, 11, 12, 13],".*FrRepCalLine.*\\[\\[(.*)\\]\\].*" ),
    ('es', 'Calendario_republicano_franc%C3%A9s', [9, 13, 17, 21], "#.*\\[\\[(.*)\\]\\]")
    ]

pattern = re.compile(".*FrRepCalLine.*\\[\\[(.*)\\]\\].*")

for (lang, title, sections, pattern_str) in titles:
    pattern = re.compile(pattern_str)
    saints = []
    for s in sections:
        # Wiki URL
        wikiUrl = f'https://{lang}.wikipedia.org/w/api.php?page={title}&action=parse&section={s}&prop=wikitext&format=json'

        rsp = requests.get(wikiUrl)
        rsp.raise_for_status()
        data = rsp.json()
        content = data['parse']['wikitext']['*']
        for line in content.split('\n'):
            match = pattern.match(line)
            if match:
                found = match.group(1)
                if '|' in found:
                    target, base = tuple(found.split('|'))
                else:
                    base = found
                    target = found
                target = target.replace(' ', '_')
                saints.append((base, target))
            elif line.startswith('#'):
                print(line)

    print(len(saints))

    with open(f'french-republican-calendar_{lang}.json', 'w') as f:
        json.dump(saints, f)

