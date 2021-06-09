# hse21_H3K27ac_ZDNA_human
# Анализ пиков гистоновой метки
## Первичная обработка файлов
- Загрузим 2 .bed файла ChIP-seq экспериментов из ENCODE (пики гистоновой метки):
wget https://www.encodeproject.org/files/ENCFF389RXK/@@download/ENCFF389RXK.bed.gz
wget https://www.encodeproject.org/files/ENCFF926NKP/@@download/ENCFF926NKP.bed.gz
- Распакуем скачанные архивы, оставив необходимые первые 5 столбцов .bed файлов:
    zcat ENCFF832EOL.bed.gz  |  cut -f1-5 > H3K4me3_A549.ENCFF832EOL.hg38.bed (аналогично для втрого файла)
- Уровень списка 1. Пункт 3.
