

import urllib.request
import re
import json

titles = [
    ('fr', 'Calendrier_r%C3%A9publicain', [29, 30, 31, 32]),
    ('en', 'French_Republican_calendar', [10, 11, 12, 13]),
    ('es', 'Calendario_republicano_francés', [9, 10, 11, 12]),
    ()
    ]

pattern = re.compile(".*FrRepCalLine.*\\[\\[(.*)\\]\\].*")

for (lang, title, sections) in titles:
    saints = []
    for s in sections:
        # Wiki URL
        wikiUrl = f'https://{lang}.wikipedia.org/w/index.php?title={title}&action=edit&section={s}'

        rsp = urllib.request.urlopen(wikiUrl)
        if not rsp.getcode() == 200:
            # Warning, Not found
            break
        for line in rsp.readlines():
            match = pattern.match(line.decode("utf-8"))
            if match:
                found = match.group(1)
                if '|' in found:
                    target, base = tuple(found.split('|'))
                else:
                    base = found
                    target = found
                target = target.replace(' ', '_')
                saints.append((base, target))

    print(saints)

    with open(f'french-republican-calendar_{lang}.json', 'w') as f:
        json.dump(saints, f)
