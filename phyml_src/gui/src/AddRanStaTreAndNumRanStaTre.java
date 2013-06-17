
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;

/**
 * Class for instantiating all the components to set the parameters whether a
 * random starting tree is added or not and the number of random starting trees.
 * 
 * @author Christoph Knapp
 * 
 */

public class AddRanStaTreAndNumRanStaTre extends JPanel implements
		ActionListener {
	/**
	 * deafult id
	 */
	private static final long serialVersionUID = 1L;
	private JComboBox addRanStaTreBox;
	private CustomTextField numRanStaTreField;
	private JLabel lab1;
	private JLabel lab2;

	/**
	 * Constructor method to instantiate all components and set their size and
	 * location.
	 * 
	 * @param isNNI
	 *            boolean : true if NNI is used for tree topology search false
	 *            otherwise
	 */
	public AddRanStaTreAndNumRanStaTre(boolean isNNI) {
		addRanStaTreBox = new JComboBox(new String[] { "no", "yes" });
		numRanStaTreField = new CustomTextField("5");
		lab1 = new JLabel("Add Random Starting Tree");
		lab2 = new JLabel("# Random Starting Trees");
		numRanStaTreField.setEnabled(false);
		if (isNNI) {
			addRanStaTreBox.setVisible(false);
			numRanStaTreField.setVisible(false);
			lab1.setVisible(false);
			lab2.setVisible(false);
		}
		addRanStaTreBox.addActionListener(this);
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		JPanel p1 = new JPanel();
		layout.setDimensions(0.44, 0.9);
		add(p1);
		layout.setDimensions(0.1, 0.9);
		add(new JPanel());
		JPanel p2 = new JPanel();
		layout.setDimensions(0.44, 0.9);
		add(p2);
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		CustomGridLayout lO1 = new CustomGridLayout();
		p1.setLayout(lO1);
		lO1.setDimensions(0.62, 1);
		p1.add(lab1);
		lO1.setDimensions(0.1, 1);
		p1.add(new JPanel());
		lO1.setDimensions(0.28, 1);
		p1.add(addRanStaTreBox);
		CustomGridLayout lO2 = new CustomGridLayout();
		p2.setLayout(lO2);
		lO2.setDimensions(0.61, 1);
		p2.add(lab2);
		lO2.setDimensions(0.1, 1);
		p2.add(new JPanel());
		lO2.setDimensions(0.23, 1);
		p2.add(numRanStaTreField);
	}

	public void actionPerformed(ActionEvent e) {
		if (addRanStaTreBox.getSelectedItem().toString().equals("no")) {
			numRanStaTreField.setEnabled(false);
		} else {
			numRanStaTreField.setEnabled(true);
		}
	}

	/**
	 * Sets the components visible if anything else but NNI is selected for Tree
	 * topology search. Otherwise all components in this panel are not visible.
	 * 
	 * @param b
	 * boolean : true if "NNI" is selected and the components are supposed to 
	 * be visible, false otherwise
	 */
	public void setIsNNI(boolean b) {
		if (b) {
			addRanStaTreBox.setVisible(false);
			numRanStaTreField.setVisible(false);
			lab1.setVisible(false);
			lab2.setVisible(false);
		} else {
			addRanStaTreBox.setVisible(true);
			numRanStaTreField.setVisible(true);
			lab1.setVisible(true);
			lab2.setVisible(true);
		}
	}

	/**
	 * Class for implementing a JTextField for Integers only.
	 * 
	 * @author Christoph Knapp
	 * 
	 */
	class CustomTextField extends JTextField implements FocusListener {
		/**
		 * default id
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * Constructor method for super() call and testing if the input String
		 * is a valid Integer input.
		 * 
		 * @param string
		 *            String : for super call to set starting input of text
		 *            field.
		 */
		public CustomTextField(String string) {
			super(string);
			if (!isNumber(string)) {
				this.setText("");
			}
			this.addFocusListener(this);
		}

		@Override
		protected Document createDefaultModel() {
			return new doubleOnlyDocument();
		}

		/**
		 * Tests wheter the parameter str is convertible in to an int value.
		 * 
		 * @param str
		 *            String : char sequence to convert.
		 * @return true if valid, false otherwise.
		 */
		private boolean isNumber(String str) {
			try {
				Integer.parseInt(str);
				return true;
			} catch (NumberFormatException e) {
				return false;
			}
		}

		/**
		 * Class extending PlainDocument to implement the document for the
		 * JTextField.
		 * 
		 * @author Christoph Knapp
		 * 
		 */
		class doubleOnlyDocument extends PlainDocument {
			/**
			 * default id
			 */
			private static final long serialVersionUID = 1L;

			@Override
			public void insertString(int offs, String str, AttributeSet a)
					throws BadLocationException {
				if (str == null || !isNumber(str)) {
					return;
				}
				super.insertString(offs, str, a);
			}
		}

		public void focusGained(FocusEvent e) {
		}

		public void focusLost(FocusEvent e) {
			if (this.getText().equals("")) {
				this.setText("5");
			}
		}
	}

	/**
	 * Whether random starting trees are used or not.
	 * 
	 * @return String : "yes" if random starting trees are used, "" otherwise.
	 */
	public String getRandStart() {
		if (addRanStaTreBox.getSelectedItem().toString().equals("yes")) {
			return "yes";
		}
		return "";
	}

	/**
	 * Retrieves the number of random start trees as a String
	 * 
	 * @return String : Either "" if no random start trees are needed or the
	 *         number as typed in by the user.
	 */
	public String getNumRandStart() {
		if (getRandStart().equals("yes")) {
			return numRanStaTreField.getText();
		}
		return "";
	}
}
