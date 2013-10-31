$( function() {
  
  // multiselect
  $("#cytokines").multiselect({
    noneSelectedText: "Cytokine combinations must contain...", 
    header: "Select Cytokines",
    selectedText: "# of # cytokines selected",
    show: ['slide', 200]  
  });
  
  // disable body scrolling when inside the multiselect
  $(".ui-multiselect-menu").mouseover( function() {
    $("body").css("overflow", "hidden");
  });
  
  $(".ui-multiselect-menu").mouseout( function() {
    $("body").css("overflow", "initial");
  });
  
  // similarly for #stats
  /*
  $("#stats").mouseover( function() {
    $("body").css("overflow", "hidden");
  });
  
  $("#stats").mouseout( function() {
    $("body").css("overflow", "initial");
  });
  */
  
  
  
});
