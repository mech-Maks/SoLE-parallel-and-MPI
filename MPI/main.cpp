#include "enviroment.h"

int main(int argc, char *argv[])
{
    /* run prog:  mpirun -n proc_amount ./a.out matrix_size block_size max_print formula_type [, file] */
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    int k;
    int p;

    MPI_Init(&argc, &argv);

    /* Разбили на процессы */
    MPI_Status status;
    MPI_Comm base_comm = MPI_COMM_WORLD;
    MPI_Comm comm;
    MPI_Group oldGroup;
    MPI_Group newGroup;
    MPI_Status st;
    MPI_Comm_size(base_comm, &p);
    MPI_Comm_rank(base_comm, &k);

    int tag = 0;
    int n = 0;
    int m = 0;
    int max_print;
    int formula_type;
    int i, j, s;
    int blocks_amount;
    int owner;
    int last = 0;
    int main_index = 0;
    int main_found = 0;
    int start_point = 0;
    int l;
    int max_blocks;
    int temp_int;
    int main_err;
    int *ranks;

    // double *buf;
    // double *a;
    // double *b;
    // double *x;

    // double *V1;
    // double *V2;
    // double *V3;

    // double *main_info;

    double min_norm = 1e50;
    double temp = 0;
    double a_norm;
    double time;

    if (argc < 5)
    {
        printf("Not enough arguments\n");
        MPI_Finalize();
        return 0;
    }

    if (argc > 6)
    {
        printf("Too much arguments\n");
        MPI_Finalize();
        return 0;
    }

    if (sscanf(argv[1], "%d", &n) != 1 ||
        sscanf(argv[2], "%d", &m) != 1 ||
        sscanf(argv[3], "%d", &max_print) != 1 ||
        sscanf(argv[4], "%d", &formula_type) != 1)
    {
        printf("Put numbers, not strings\n");
        MPI_Finalize();
        return 0;
    }

    if (formula_type == 0 && argc < 6)
    {
        printf("No input-file\n");
        MPI_Finalize();
        return 0;
    }

    if (n <= 0 ||
        m <= 0 ||
        max_print < 0 ||
        formula_type < 0)
    {
        printf("Put in possitive numbers\n");
        MPI_Finalize();
        return 0;
    }

    if (formula_type < 0 || formula_type > 4)
    {
        printf("Unknown formula type\n");
        MPI_Finalize();
        return 0;
    }

    if (n < m)
        m = n;
    if (n < max_print)
        max_print = n;

    last = (n % m == 0 ? 0 : 1);

    if (n / m + last < p)
    {
        temp_int = p;
        while (n / m + last < temp_int)
        {
            temp_int--;
        }

        ranks = new int[temp_int];
        for (i = 0; i < temp_int; i++)
            ranks[i] = i;

        MPI_Comm_group(base_comm, &oldGroup);
        MPI_Group_incl(oldGroup, temp_int, ranks, &newGroup);
        MPI_Comm_create(base_comm, newGroup, &comm);
        p = temp_int;

        delete[] ranks;
        if (k >= temp_int)
        { // закрываем ненужные процессы
            MPI_Finalize();
            return 0;
        }
    }
    else
    {
        comm = MPI_COMM_WORLD;
    }

    double *a = new double[n * m * get_max_rows_block(n, m, p)];
    double *b = new double[1 * m * get_max_rows_block(n, m, p)];
    int *x = new int[n];
    double *buf = new double[n * m];
    double *b1 = new double[m];
    double *b2 = new double[m];
    double *b3 = new double[m];
    double *b_buf = new double[m];

    for (i = 0; i < 1 * m * get_rows_block(n, m, p, k); i++)
        b[i] = 0;

    for (i = 0; i < n * m; i++)
        buf[i] = 0;

    if (!a)
    {
        printf("Computer is dead\n");
        MPI_Finalize();
        return 0;
    }

    // занулили матрицу
    for (i = 0; i < n * m * get_max_rows_block(n, m, p); i++)
        a[i] = 0;

    for (i = 0; i < n; i++)
        x[i] = i;

    // V1 = (double *)calloc(m * m, sizeof(double));
    // V2 = (double *)calloc(m * m, sizeof(double));
    // V3 = (double *)calloc(m * m, sizeof(double));

    double *V1 = new double[m * m];
    double *V2 = new double[m * m];
    double *V3 = new double[m * m];
    for (i = 0; i < m * m; i++)
    {
        V1[i] = 0;
        V2[i] = 0;
        V3[i] = 0;
    }

    // main_info = (double *)calloc(3, sizeof(double));
    double *main_info = new double[3];
    for (i = 0; i < 3; i++)
        main_info[i] = 0;

    if (formula_type != 0)
    {
        init_matrix(a, n, m, p, k, formula_type, &formula);
    }
    else
    {
        if ((read_matrix(a, n, m, p, k, argv[5], buf) != 0))
        {
            printf("Failed while reading matrix from file\n");

            delete[] V1;
            delete[] V2;
            delete[] V3;
            delete[] a;
            delete[] buf;
            delete[] x;
            delete[] b;
            delete[] main_info;
            delete[] b1;
            delete[] b2;
            delete[] b3;
            delete[] b_buf;

            MPI_Finalize();
            return 0;
        }
    }

    // Заполняем b-вектор
    for (i = 0; i < get_rows_block(n, m, p, k) * m; i++)
    {
        b[i] = 0;

        for (j = 0; j < n; j++)
        {
            if (j % 2 == 0)
            {
                b[i] += a[i * n + j];
            }
            else
            {
                continue;
            }
        }
    }

    print_matrix(a, n, m, p, k, buf, max_print);

    /*
    * Заполнение, чтение и печать корректны, приступаем к вычислениям
    */

    time = MPI_Wtime();

    blocks_amount = n / m;
    a_norm = matrixNorm(a, m * get_rows_block(n, m, p, k), n);

    if (last)
    {
        l = n - (n / m) * m;
    }

    // if (k == 1)
    //     print_fixed_array(b, get_rows_block(n, m, p, k) * m, 1);

    MPI_Barrier(comm);

    /*================================== ПРЯМОЙ ХОД ==================================*/

    for (s = 0; s < blocks_amount; s++)
    { // шаг прямого хода
        owner = s % p;
        main_err = 1; // 1 - в строке не нашелся обратимый блок, 0 - нашелся
        min_norm = 1e50;
        main_index = 0;
        main_found = 0;

        if (k == owner)
        { // находимся в процессе, владеющим данной блочной строкой
            memcpy(buf, a + (s / p) * m * n, n * m * 8);
        }

        MPI_Bcast(buf, n * m, MPI_DOUBLE, owner, comm);

        i = k;
        while (i < s)
        {
            i += p;
        }

        while (i < blocks_amount)
        { // Каждый процесс ищет главный блок в своих кусках выделенных внутри буфферов
            temp = RecipMatrixNorm(buf, n, m, i, a_norm, V1, V3);
            if (temp > 0)
            {
                if (min_norm > temp)
                {
                    min_norm = temp;
                    main_index = i;
                    // printf("Чет нашли в процессе %d\n", k);
                }

                main_found = 1;
                main_err = 0;
            }
            i += p;
        }

        main_info[0] = (double)main_found;
        main_info[1] = min_norm;
        main_info[2] = (double)main_index;

        MPI_Barrier(comm);

        if (k != owner)
        { // Отправляем информацию из других процессов
            MPI_Send(main_info, 3, MPI_DOUBLE, owner, tag, comm);
        }
        else
        { // Принимаем информацию из других процессов
            for (i = 0; i < p; i++)
            {
                if (i == owner)
                {
                    continue;
                }

                MPI_Recv(main_info, 3, MPI_DOUBLE, i, tag, comm, &st);

                if ((int)main_info[0] > 0)
                { // В i процессе нашелся хотя бы 1 обратимый блок
                    temp = main_info[1];

                    if (temp < min_norm)
                    {
                        main_index = (int)main_info[2];
                    }

                    main_err = 0;
                }
            }
        }

        MPI_Barrier(comm);

        MPI_Bcast(&main_index, 1, MPI_INT, owner, comm);
        MPI_Bcast(&main_err, 1, MPI_INT, owner, comm);

        if (main_err)
        { // В строке не нашлось обратимого блока
            if (k == 0)
            {
                printf("на %d шаге не нашлось обратимого блока в строке\n", s);
            }

            delete[] V1;
            delete[] V2;
            delete[] V3;
            delete[] a;
            delete[] buf;
            delete[] x;
            delete[] b;
            delete[] main_info;
            delete[] b1;
            delete[] b2;
            delete[] b3;
            delete[] b_buf;

            MPI_Finalize();
            return 0;
        }

        if (main_index > s && k == 0)
        { // Меняем блоки в строках
            // printf("Меняем %d и %d!\n", s, main_index);
            for (i = 0; i < m; i++)
            {
                temp_int = x[s * m + i];
                x[s * m + i] = x[main_index * m + i];
                x[main_index * m + i] = temp_int;
            }
        }

        // В каждом процессе надо поменять местами s и main_index столбцы блоков
        // В главном не забыть поменять индексы у вектора X

        if (main_index > s)
            swap_columns(a, buf, n, m, p, k, s, main_index);

        copyBlockString(buf, V1, s, n, m, m, m);
        recipMatrix(V1, V2, m, a_norm); /* обратная лежит в V2 */

        setIdentityMatrix(V1, m, m);

        if (k == owner)
        { // домножаем в процессе, который владеет вычитаемой строкой
            writeBlockMatix(a, V1, start_point, s, n, m, m, m);
            writeBlockString(buf, V1, s, n, m, m, m);
            for (i = s + 1; i < n / m; i++)
            {
                copyBlockMatix(a, V1, start_point, i, n, m, m, m);
                multiplyMatricesQuadratic(V2, V1, V3, m);
                // multiplyMatrices(V2, m, m, V1, m, m, V3);

                writeBlockMatix(a, V3, start_point, i, n, m, m, m);
                writeBlockString(buf, V3, i, n, m, m, m);
            }

            if (last)
            {
                copyBlockMatix(a, V1, start_point, n / m, n, m, m, l);
                multiplyMatrices(V2, m, m, V1, m, l, V3);

                writeBlockMatix(a, V3, start_point, n / m, n, m, m, l);
                writeBlockString(buf, V3, n / m, n, m, m, l);
            }

            copyBlockVector(b, b1, start_point, m, m);
            multiplyMatrices(V2, m, m, b1, m, 1, b2);
            writeBlockVector(b, b2, start_point, m, m);

            copyBlockVector(b, b_buf, start_point, m, m);

            start_point += 1;
        }
        else
        {
            setIdentityMatrix(V1, m, m);
            writeBlockString(buf, V1, s, n, m, m, m);
            for (i = s + 1; i < n / m; i++)
            {
                copyBlockString(buf, V1, i, n, m, m, m);
                multiplyMatricesQuadratic(V2, V1, V3, m);
                // multiplyMatrices(V2, m, m, V1, m, m, V3);

                writeBlockString(buf, V3, i, n, m, m, m);
            }

            if (last)
            {
                copyBlockString(buf, V1, n / m, n, m, m, l);
                multiplyMatrices(V2, m, m, V1, m, l, V3);

                writeBlockString(buf, V3, n / m, n, m, m, l);
            }
        }

        MPI_Barrier(comm);
        MPI_Bcast(b_buf, m, MPI_DOUBLE, owner, comm);

        // Вычитаем строки
        max_blocks = get_rows_block(n, m, p, k);
        for (i = start_point; i < max_blocks; i++)
        {
            copyBlockMatix(a, V1, i, s, n, m, m, m);

            for (j = s; j < n / m; j++)
            {
                copyBlockString(buf, V2, j, n, m, m, m);
                multiplyMatricesQuadratic(V1, V2, V3, m);
                // multiplyMatrices(V1, m, m, V2, m, m, V3);

                copyBlockMatix(a, V2, i, j, n, m, m, m);
                subtractMatrices(V2, V3, m, m);

                writeBlockMatix(a, V2, i, j, n, m, m, m);
            }

            if (last)
            {
                copyBlockString(buf, V2, n / m, n, m, m, l);
                multiplyMatrices(V1, m, m, V2, m, l, V3);

                copyBlockMatix(a, V2, i, n / m, n, m, m, l);
                subtractMatrices(V2, V3, m, l);

                writeBlockMatix(a, V2, i, j, n, m, m, l);
            }

            copyBlockVector(b, b1, i, m, m);
            multiplyMatrices(V1, m, m, b_buf, m, 1, b2);
            subtractMatrices(b1, b2, m, 1);
            writeBlockVector(b, b1, i, m, m);
        }

        MPI_Barrier(comm);
    }

    // Обработка остаточного элемента в правом нижнем углу

    MPI_Barrier(comm);

    if (last > 0 && k == ((n / m) % p))
    {
        i--;
        main_err = 1;
        copyBlockMatix(a, V1, n / m / p, n / m, n, m, l, l);
        temp = recipMatrix(V1, V2, l, 10 * a_norm);

        if (temp > 0)
        {
            main_err = 0;
            copyBlockVector(b, b1, i, m, l);

            multiplyMatrices(V2, l, l, b1, l, 1, b2);
            writeBlockVector(b, b2, i, m, l);

            setIdentityMatrix(V1, l, l);
            writeBlockMatix(a, V1, n / m / p, n / m, n, m, l, l);
        }
    }

    MPI_Bcast(&main_err, 1, MPI_INT, (n / m) % p, comm);
    if (main_err)
    { // Элемент в правом нижнем углу необратим
        // printf("ОШИБКА\n");
        if (k == 0)
        {
            printf("Элемент в правом нижнем углу необратим\n");
        }

        delete[] V1;
        delete[] V2;
        delete[] V3;
        delete[] a;
        delete[] buf;
        delete[] x;
        delete[] b;
        delete[] main_info;
        delete[] b1;
        delete[] b2;
        delete[] b3;
        delete[] b_buf;

        MPI_Finalize();
        return 0;
    }

    MPI_Barrier(comm);

    // прямой ход вроде правильный

    /* ================================== ОБРАТНЫЙ ХОД ================================== */

    start_point = get_rows_block(n, m, p, k) - 1;

    // обработка последнего столбца
    if (last > 0)
    {
        if (k == (n / m) % p)
        { // Процесс, владеющий последней строкой
            copyBlockVector(b, b1, n / m / p, m, l);
            start_point -= 1;
        }

        MPI_Bcast(b1, l, MPI_DOUBLE, (n / m) % p, comm);

        for (i = start_point; i >= 0; i--)
        {
            copyBlockMatix(a, V1, i, n / m, n, m, m, l);

            multiplyMatrices(V1, m, l, b1, l, 1, b2);

            copyBlockVector(b, b3, i, m, m);
            // print_fixed_array(b3, l, 1);

            subtractMatrices(b3, b2, m, 1);

            writeBlockVector(b, b3, i, m, m);

            setZeroMatrix(V1, m, l);
            writeBlockMatix(a, V1, i, n / m, n, m, m, l);
        }
    }

    for (s = (n / m - 1); s >= 0; s--)
    { // Шаг обратного хода
        owner = s % p;

        if (k == owner)
        { // Раздаем блок вектора b
            copyBlockVector(b, b1, s / p, m, m);
            start_point -= 1;
        }

        MPI_Bcast(b1, m, MPI_DOUBLE, owner, comm);

        for (j = start_point; j >= 0; j--)
        {
            copyBlockMatix(a, V1, j, s, n, m, m, m);
            multiplyMatrices(V1, m, m, b1, m, 1, b2);

            copyBlockVector(b, b3, j, m, m);
            subtractMatrices(b3, b2, m, 1);

            writeBlockVector(b, b3, j, m, m);

            setZeroMatrix(V1, m, m);
            writeBlockMatix(a, V1, j, s, n, m, m, m);
        }
    }

    MPI_Barrier(comm);
    print_matrix(a, n, m, p, k, buf, max_print);

    // for (int h = 0; h < p; h++)
    // {
    //     if (k == h)
    //         print_fixed_array(b, m * get_rows_block(n, m, p, k), 1);

    //     MPI_Barrier(comm);
    // }

    time = MPI_Wtime() - time;

    // ================================== ВРЕМЯ СЧИТАТЬ НЕВЯЗКУ ==================================

    /* 
     * Для невязки каждому процессу требуется буффер под полученные значения и
     * под вектор b. 
     * Собираем полученный вектор b в нулевом процессе (главном) и формируем
     * ответ (вектор values)
     * 
     */

    /* Собираем все значения в главном процессе */
    double *values = new double[1 * m * (n / m + last)];
    double *b_total = new double[1 * m * (n / m + last)];

    if (k == 0)
    {
        for (i = 0; i < get_rows_block(n, m, p, k); i++)
        {
            for (j = 0; j < m; j++)
            {
                b_total[i * m * p + k * m + j] = b[i * m + j];
            }
        }

        for (i = 1; i < p; i++)
        {
            MPI_Recv(b, m * get_rows_block(n, m, p, i), MPI_DOUBLE, i, tag, comm, &st);

            for (int h = 0; h < get_rows_block(n, m, p, i); h++)
            {
                for (j = 0; j < m; j++)
                {
                    b_total[h * m * p + i * m + j] = b[h * m + j];
                }
            }
        }
    }
    else
    {
        MPI_Send(b, get_rows_block(n, m, p, k) * m, MPI_DOUBLE, 0, tag, comm);
    }

    MPI_Barrier(comm);

    // Получили полный вектор с.з.
    if (k == 0)
    {
        // printf("Total b:\n");
        // print_fixed_array(b_total, m * (n / m + last), 1);

        for (i = 0; i < n; i++)
        {
            values[x[i]] = b_total[i];
        }

        // print_fixed_array(values, m * (n / m + last), 1);
    }

    // Восстановим матрицу а и векторы b

    if (formula_type != 0)
    {
        init_matrix(a, n, m, p, k, formula_type, &formula);
    }
    else
    {
        read_matrix(a, n, m, p, k, argv[5], buf);
    }

    // Заполняем b-вектор
    for (i = 0; i < get_rows_block(n, m, p, k) * m; i++)
    {
        b[i] = 0;

        for (j = 0; j < n; j++)
        {
            if (j % 2 == 0)
            {
                b[i] += a[i * n + j];
            }
            else
            {
                continue;
            }
        }
    }

    // Формируем полностью вектор b в главном процессе (затем раздадим его)
    if (k == 0)
    {
        for (i = 0; i < get_rows_block(n, m, p, k); i++)
        {
            for (j = 0; j < m; j++)
            {
                b_total[i * m * p + k * m + j] = b[i * m + j];
            }
        }

        for (i = 1; i < p; i++)
        {
            MPI_Recv(b, m * get_rows_block(n, m, p, i), MPI_DOUBLE, i, tag, comm, &st);

            for (int h = 0; h < get_rows_block(n, m, p, i); h++)
            {
                for (j = 0; j < m; j++)
                {
                    b_total[h * m * p + i * m + j] = b[h * m + j];
                }
            }
        }

        // print_fixed_array(b_total, m * (n / m + last), 1);
    }
    else
    {
        MPI_Send(b, get_rows_block(n, m, p, k) * m, MPI_DOUBLE, 0, tag, comm);
    }

    // Раздаем вектор b_total и values другим процессам, чтобы те использовали его при подсчете невязки
    MPI_Bcast(b_total, m * (n / m + last), MPI_DOUBLE, 0, comm);
    MPI_Bcast(values, m * (n / m + last), MPI_DOUBLE, 0, comm);

    // if (k == 0)
    //     print_fixed_array(values, m * (n / m + last), 1);

    // Начинаем перемножать и заносить результат в массив b
    for (i = 0; i < get_rows_block(n, m, p, k); i++)
    {
        setZeroMatrix(b2, m, 1);

        for (j = 0; j < n / m; j++)
        {
            copyBlockMatix(a, V1, i, j, n, m, m, m);
            copyBlockVector(values, b_buf, j, m, m);

            multiplyMatrices(V1, m, m, b_buf, m, 1, b1);
            addMatrices(b2, b1, m, 1);
        }

        if (last)
        {
            copyBlockMatix(a, V1, i, n / m, n, m, m, l);
            copyBlockVector(values, b_buf, n / m, m, l);

            multiplyMatrices(V1, m, l, b_buf, l, 1, b1);
            addMatrices(b2, b1, m, 1);
        }

        writeBlockVector(b, b2, i, m, m);
    }

    // Собираем промежуточные умножения в главном процессе
    if (k == 0)
    {
        for (i = 0; i < get_rows_block(n, m, p, k); i++)
        {
            for (j = 0; j < m; j++)
            {
                values[i * m * p + k * m + j] = b[i * m + j];
            }
        }

        for (i = 1; i < p; i++)
        {
            MPI_Recv(b, m * get_rows_block(n, m, p, i), MPI_DOUBLE, i, tag, comm, &st);

            for (int h = 0; h < get_rows_block(n, m, p, i); h++)
            {
                for (j = 0; j < m; j++)
                {
                    values[h * m * p + i * m + j] = b[h * m + j];
                }
            }
        }
    }
    else
    {
        MPI_Send(b, get_rows_block(n, m, p, k) * m, MPI_DOUBLE, 0, tag, comm);
    }

    if (k == 0)
    {
        subtractMatrices(values, b_total, n, 1);
        temp = matrixNorm(values, n, 1) / matrixNorm(b_total, n, 1);

        printf("%s : residual = %e elapsed = %.2f for s = %d n = %d m = %d p = %d\n\n",
               argv[0], temp, time, formula_type, n, m, p);
    }

    delete[] V1;
    delete[] V2;
    delete[] V3;
    delete[] a;
    delete[] buf;
    delete[] x;
    delete[] b;
    delete[] main_info;
    delete[] b1;
    delete[] b2;
    delete[] b3;
    delete[] b_buf;
    delete[] values;
    delete[] b_total;

    (void)status;

    MPI_Finalize();
    return 0;
}