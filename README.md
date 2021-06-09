# hse21_H3K27ac_ZDNA_human
# Анализ пиков гистоновой метки
## Первичная обработка файлов
- Загрузим 2 .bed файла ChIP-seq экспериментов из ENCODE (пики гистоновой метки):
 ```
 wget https://www.encodeproject.org/files/ENCFF389RXK/@@download/ENCFF389RXK.bed.gz
 wget https://www.encodeproject.org/files/ENCFF926NKP/@@download/ENCFF926NKP.bed.gz
 ```
- Распакуем скачанные архивы, оставив необходимые первые 5 столбцов .bed файлов:
```
zcat ENCFF832EOL.bed.gz  |  cut -f1-5 > H3K4me3_A549.ENCFF832EOL.hg38.bed (аналогично для второго файла)
```
- Переведем координаты ChIP-seq пиков к версии генома hg19 с помощью liftOver:
     - Загрузим необходимый файл для конверации hg38 в hg19:
    ```
    wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz
    ```
    - Запустим liftOver на файлах .bed:
    >Используем команду в формате: *liftOver input.bed hg18ToHg19.over.chain.gz output.bed unlifted.bed*,
    >>где input.bed - файл с исходными данными, hg18ToHg19.over.chain.gz - вспомогательный файл для конвертации,
    >>output.bed - файл, содержащий конвертированные координаты, unlifted.bed будет содержать все позиции генома, которые нельзя перевести в нужную вурсию генома
     ```
    liftOver   H3K4me3_A549.ENCFF832EOL.hg38.bed   hg38ToHg19.over.chain.gz   H3K4me3_A549.ENCFF832EOL.hg19.bed   H3K4me3_A549.ENCFF832EOL.unmapped.bed (аналогично для второго     файла)
    ```
- Построим гистограмму длин участков для каждого эксперимента до и после конвертации к нужной версии генома:
  > Для обработки данных используем скрипт */src/filter_peaks.R*
    - ***Гистограммы для версии генома hg38***
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF389RXK.hg38.init.hist.png) 
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF926NKP.hg38.init.hist.png)
    - ***Гистограммы для версии генома hg19***
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF389RXK.hg19.init.hist.png) 
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF926NKP.hg19.init.hist.png)
- Проведем фильтрацию пиков по длине:
 >Установим порог 5kb (не длиннее 5000 пар оснований)
 >>Для обработки данных используем скрипт */src/filter_peaks.R - раздел Remove long peaks*
    - ***Гистограммы для filter_peaks***
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF389RXK.hg19.filtered.hist.png) 
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF926NKP.hg19.filtered.hist.png)
- Статистика по количеству пиков:
  |  | ENCFF926NKP | ENCFF389RXK |
  |----------------|:---------:|:---------:|
  | Для версии генома ***hg38*** | 118990 | 116220 |
  | Для версии генома ***hg19*** | 118767 | 115997 | 
  | Для ***filter_peaks*** | 115115 | 112077 |  

    
