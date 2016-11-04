$(document).ready(function() {

    $(".hideshow").each(function(){
        var theid=$(this).attr('id');
        $("#" + theid).click(function(){
            $(this).parent().css('width', $(this).width());
            $("#table-" + theid).toggle();
        });
    });
    // Go through each table.obstable object
    $("table.obstable").each(function(){

        // get the id of the current table.
        var theid = $(this).attr('id');
        console.log("the id is " + theid);

        // replace the tr#filterrow th in the header with search inputs
        $('#' + theid + ' thead tr#filterrow th').each( function () {
            var title = $(this).text();
            $(this).html( '<input type="text" placeholder="Search '+title+'" />' );
        } );

        // Initialise the table
        var table=$('#'+theid).DataTable( {
            dom: 'Biprtip',
            buttons: [
                'colvis',
                'csvHtml5',
            ],
             "columnDefs": [
                 {targets: ["instrument","trans"], visible: false }],
            colReorder: true,
            "paging": false,
            orderCellsTop: true,
            fixedHeader: true,
            orderClasses: false,
        } )
           // Apply the search
        $("table#"+theid+" thead input").on( 'keyup change', function() {
            table
                .column( $(this).parent().index()+':visible' )
                .search( this.value )
                .draw();
        } );

    });
    // // DataTable
    // var table = $('table.display').DataTable( {
    //     dom: 'Bipfrtip',
    //     buttons: [
    //         'colvis',
    //         'csvHtml5',
    //     ],
    //     colReorder: true,
    //     "paging": false,
    //     orderCellsTop: true,
    //     fixedHeader: true,
    // } );



});
