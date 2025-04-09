import os
import re
from datetime import datetime


def print_colored_ascii(html_path):
    try:
        with open(html_path, 'r', encoding='utf-8') as file:
            content = file.read()

            pattern = r'<div style="margin: 20px 0;">(.*?)</div>'
            match = re.search(pattern, content, re.DOTALL)

            if match:
                art_content = match.group(1)
                art_content = art_content.replace('<br/>', '\n').replace(
                    '<br>', '\n')

                art_content = re.sub(
                    r'<span style="color: rgb\((\d+),\s*(\d+),\s*(\d+)\)">(.*?)</span>',
                    lambda
                        m: f"\033[38;2;{m.group(1)};{m.group(2)};{m.group(3)}m{m.group(4)}\033[0m",
                    art_content
                )

                art_content = re.sub(r'<[^>]+>', '', art_content)
                art_content = art_content.replace('&nbsp;', ' ')
                lines = art_content.split('\n')
                formatted_lines = []
                for line in lines:
                    if line.strip():
                        if not line.endswith('\033[0m'):
                            line += '\033[0m'
                        formatted_lines.append(line)

                print('\n\n')
                print('\n'.join(formatted_lines))
                print('\033[0m')

    except Exception as e:
        print(f"Error al leer el archivo HTML: {e}")



if __name__ == "__main__":
    current_dir = os.path.dirname(os.path.abspath(__file__))
    html_path = os.path.join(current_dir, "binary_imap.html")
    print_colored_ascii(html_path)
