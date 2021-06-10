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
## Анализ длины пиков
- Построим гистограмму длин участков для каждого эксперимента до и после конвертации к нужной версии генома:
  > Для обработки данных используем скрипт */src/filter_peaks.R*
    - ***Гистограммы для версии генома hg38***
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF389RXK.hg38.init.hist.png) 
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF926NKP.hg38.init.hist.png)
    - ***Гистограммы для версии генома hg19***
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF389RXK.hg19.init.hist.png) 
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/filter_peaks.H3K27ac_A549.ENCFF926NKP.hg19.init.hist.png)
- Проведем фильтрацию пиков по длине:
    - Рассмотрим и сравним самые длинные пики для разных версий геномов:
      >|  | Версия генома hg38 |  Версия генома hg19 |
      >|----------------|:---------------------------------------------:|:------------------------------------------------:| 
      >| ENCFF926NKP|   chrom    start       end      name score   len  |   chrom     start       end      name score   len |
      >|            | 1 chr20  47299098  47363843  Peak_191  1000 64745 | 1 chr20  45927842  45992587  Peak_191  1000 64745 |
      >|            | 2 chr16  87821385  87872899 Peak_4576  1000 51514 | 2 chr16  87854991  87906505 Peak_4576  1000 51514 |
      >|            | 3 chr3  53226233  53275734 Peak_1352  1000 49501  | 3  chr3  53260249  53309750 Peak_1352  1000 49501 |
      >|            | 4 chr15 98849239  98898727 Peak_1710  1000 49488  | 4 chr15  99392468  99441956 Peak_1710  1000 49488 |
      >|            | 5 chr22 30219891  30262948  Peak_660  1000 43057  | 5 chr22  30615880  30658937  Peak_660  1000 43057 |
      >|            | 6 chr5 172894105 172935603   Peak_86  1000 41498  | 6  chr5 172321108 172362606   Peak_86  1000 41498 |
      >| ENCFF389RXK|   chrom    start      end      name score   len   |   chrom    start      end      name  score   len  |
      >|            | 1 chr20 47301248 47363833   Peak_14  1000 62585   | 1 chr20 45929992 45992577   Peak_14  1000 62585   |
      >|            | 2 chr20 50411182 50465557  Peak_724  1000 54375   | 2 chr20 49027719 49082094  Peak_724  1000 54375   |
      >|            | 3 chr16 87819289 87872923 Peak_1120  1000 53634   | 3 chr16 87852895 87906529 Peak_1120  1000 53634   |
      >|            | 4  chr7 47567332 47617458   Peak_82  1000 50126   | 4  chr7 47606930 47657056   Peak_82  1000 50126   |
      >|            | 5 chr15 98849850 98898737  Peak_863  1000 48887   | 5 chr15 99393079 99441966  Peak_863  1000 48887   |
      >|            | 6 chr22 30219910 30264229  Peak_270  1000 44319   | 6 chr22 30615899 30660218  Peak_270  1000 44319   |
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
## Анализ расположения пиков относительно аннотированных генов
- Рассмотрим, где располагаются пики гистоновой метки относительно аннотированных генов:
    - Строим график типа пай-чарт (ENCFF389RXK - ENCFF926NKP):
      > Для обработки данных используем скрипт */src/chip_seeker.R.R*
   ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/chip_seeker.H3K27ac_A549.ENCFF389RXK.hg19.filtered.plotAnnoPie.png) 
   ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/chip_seeker.H3K27ac_A549.ENCFF926NKP.hg19.filtered.plotAnnoPie.png)
## Объеденение отфильтрованных ChIP-seq пиков с помощью утилиты bedtools merge
- Запустим команду, предварительно отсортировав единый .bed файл:
   ```
   cat  *.filtered.bed  |   sort -k1,1 -k2,2n   |   bedtools merge   >  H3K4me3_A549.merge.hg19.bed 
   ```
## Визуализация данных в геномном браузере и проверяем корректность работы bedtools merge
- Визуализируем исходные два набора ChIP-seq пиков, а также их объединение в UCSC Genome Browser:
    - На сайте UCSC Genome Browser:
      ```
      Выбираем версию генома (вкладка Genomes) => MyData => Custom Tracks
       ```
    - Добавим сustom Tracks командами:
      ### **Для filtered.bed файлов**
         ```
         track visibility=dense name="ENCFF926NKP" description="H3K27ac_A549.ENCFF926NKP.hg19.filtered.bed"
         https://raw.githubusercontent.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/main/data/H3K27ac_A549.ENCFF926NKP.hg19.filtered.bed

         track visibility=dense name="ENCFF389RXK" description="H3K27ac_A549.ENCFF389RXK.hg19.filtered.bed"
         https://raw.githubusercontent.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/main/data/H3K27ac_A549.ENCFF389RXK.hg19.filtered.bed
         ```
      ### **Для объединеного файла**
          ```
          track visibility=dense name="ChIP_merge"  color=50,50,200   description="H3K27ac_A549.merge.hg19.bed"
          https://raw.githubusercontent.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/main/data/H3K27ac_A549.merge.hg19.bed
          ```
     - Ссылка на визуализированные данные:
       ```
      http://genome.ucsc.edu/s/KamillaZyal/ENCFF926NKP_ENCFF389RXK_filtered
       ```
     - Сделаем вывод о корректности работы bedtools merge:
       ```
       Ниже представлен скриншот участка с координатами chr3:187,432,810-187,482,809, подтверждающий корректную работу утилиты bedtools merge
       ```
       ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/mergeUCSC.png) 
# Анализ участков вторичной стр-ры ДНК
## Загрузка файл со вторичной стр-рой ДНК
- Загружаем файл соответствующий предсказания ZDNA_DeepZ:
```
 wget https://github.com/Nazar1997/DeepZ/blob/master/annotation/DeepZ.bed
 ```
 ## Потроение диаграмм
 - Построим гистограмму длин участков:
     > Для обработки данных используем скрипт */src/len_hist.R*
     - ***Гистограммы для DeepZ***
       ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/len_hist.DeepZ.png) 
 - Рассмотрим, где располагаются пики гистоновой метки относительно аннотированных генов:
       > Для обработки данных используем скрипт */src/chip_seeker.R.R*
     - Строим график типа пай-чарт:
     ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/chip_seeker.DeepZ.plotAnnoPie.png)
# Анализ пересечений гистоновой метки и стр-ры ДНК
## Первичная обработка
- Находие пересечения гистоновой меткой и стр-рами ДНК:
   ```
   bedtools intersect  -a DeepZ.bed   -b  H3K27ac_A549.merge.hg19.bed  >  H3K427ac_A549.intersect_with_DeepZ.bed 
   ```
- Построим гистограмму длин участков:
  > Для обработки данных используем скрипт */src/len_hist.R*
    - ***Гистограммы для H3K427ac_A549.intersect_with_DeepZ***
    ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/len_hist.H3K27ac_A549.intersect_with_DeepZ.png) 
## Визуализация данных в геномном браузере 
- Визуализируем исходные два набора ChIP-seq пиков, а также их объединение в UCSC Genome Browser:
    - На сайте UCSC Genome Browser:
      ```
      Выбираем версию генома (вкладка Genomes) => MyData => Custom Tracks
       ```
    - Добавим сustom Tracks командами:
      ### **Для DeepZ.bed файлов**
         ```
         track visibility=dense name="DeepZ"  color=0,200,0  description="DeepZ"
         https://raw.githubusercontent.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/main/data/DeepZ.bed
         ```
      ### **Для intersect_with_DeepZ.bed файла**
          ```
          track visibility=dense name="intersect_with_DeepZ"  color=200,0,0  description="H3K27ac_A549.intersect_with_DeepZ.bed"
          https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/raw/main/data/H3K27ac_A549.intersect_with_DeepZ.bed
          ```
     - Ссылка на визуализированные данные:
       ```
       http://genome.ucsc.edu/s/KamillaZyal/hg19
       ```
     - Найдем пересечения между гистоновой меткой и стр-рой ДНК (рядом с аннотированным геном):
        ### **Координаты  chr3:105,085,585-105,085,615**
        ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/intersect1.png) 
        ### **Координаты chr14:105,877,676-105,877,817**
        ![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/intersect2.png)
## Ассоция полученных пересечений с ближайшими генами
- Используем R-библиотеку ChIPpeakAnno для ассоциаций пересееней (диапозон: -2000-3000):
  > Для обработки данных используем скрипт */src/len_hist.R*
 ```
 Получаем файлы: [H3K27ac_A549.intersect_with_DeepZ.genes_uniq](http://pantherdb.org/)
 ```
 [H3K27ac_A549.intersect_with_DeepZ.genes_uniq](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/data/H3K27ac_A549.intersect_with_DeepZ.genes_uniq.txt)
 [H3K27ac_A549.intersect_with_DeepZ.genes](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/data/H3K27ac_A549.intersect_with_DeepZ.genes.txt)
## GO-анализ
- Для этого можно воспользоваться сайтом [http://pantherdb.org/](http://pantherdb.org/):
  > Загрузим список геном из H3K27ac_A549.intersect_with_DeepZ.genes_uniq
- Наиболее значимые категории (минимальные значения FDR):
    > Полный список [pantherdb_GO_analysis](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/data/pantherdb_GO_analysis.txt.txt)
![Alt-текст](https://github.com/KamillaZyal/hse21_H3K27ac_ZDNA_human/blob/main/images/GOan.png)
    
