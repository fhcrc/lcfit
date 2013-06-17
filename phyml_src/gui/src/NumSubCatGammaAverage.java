
import java.awt.Color;
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
 * Implements all components to specify the "# Subcategories", 
 * "Gamma", "Alpha", and the "Average".
 * 
 * @author Christoph Knapp
 *
 */
public class NumSubCatGammaAverage extends JPanel implements ActionListener,
		FocusListener {
	/**
	 * default id
	 */
	private static final long serialVersionUID = 1L;
	public CustomTextField subCat;
	private JComboBox gammaBox;
	private CustomTextField alpha;
	private JComboBox average;
	private JLabel lab2;
	private JLabel lab1;
	private JLabel lab3;
	private JLabel lab4;

	/**
	 * Constructor method implements all components to specify the "# Subcategories", "Gamma", 
	 * "Alpha", and the "Average", it instantiates all components, adds listeners sets size 
	 * and position of components and adds the components.
	 */
	public NumSubCatGammaAverage() {
		subCat = new CustomTextField("4", false);
		gammaBox = new JComboBox(new String[] { "estimated", "fixed" });
		gammaBox.addActionListener(this);
		alpha = new CustomTextField("", true);
		alpha.setEnabled(false);
		alpha.addFocusListener(this);
		average = new JComboBox(new String[] { "mean", "median" });
		lab1 = new JLabel("# of rate classes");
		lab2 = new JLabel("Shape parameter");
		lab3 = new JLabel("Value");
		lab4 = new JLabel("Average");
		CustomGridLayout layout = new CustomGridLayout();
		setLayout(layout);
		layout.setDimensions(1, 0.1);
		add(new JPanel());
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		layout.setDimensions(0.25, 0.9);

		JPanel p1 = new JPanel();
		add(p1);
		//layout.setDimensions(0.01, 0.8);
		//add(new JPanel());
		layout.setDimensions(0.5, 0.9);

		JPanel p2 = new JPanel();
		//p2.setBackground(Color.blue);
		add(p2);
		//layout.setDimensions(0.01, 0.8);
		//add(new JPanel());
		layout.setDimensions(0.23, 0.9);
		JPanel p3 = new JPanel();
		add(p3);
		layout.setDimensions(0.01, 0.9);
		add(new JPanel());
		//layout.setDimensions(1, 0.1);
		//add(new JPanel());
		CustomGridLayout lO1 = new CustomGridLayout();
		p1.setLayout(lO1);
		lO1.setDimensions(0.7, 1);
		p1.add(lab1);
		lO1.setDimensions(0.25, 1);
		p1.add(subCat);
		lO1.setDimensions(0.05, 1);
		p1.add(new JPanel());
		CustomGridLayout lO2 = new CustomGridLayout();
		p2.setLayout(lO2);
		lO2.setDimensions(0.36, 1);
		p2.add(lab2);
		lO2.setDimensions(0.25, 1);
		p2.add(gammaBox);
		lO2.setDimensions(0.03, 1);
		p2.add(new JPanel());
		lO2.setDimensions(0.18, 1);
		p2.add(lab3);
		lO2.setDimensions(0.18, 1);
		p2.add(alpha);
		CustomGridLayout lO3 = new CustomGridLayout();
		p3.setLayout(lO3);
		lO3.setDimensions(0.04, 1);
		p3.add(new JPanel());
		lO3.setDimensions(0.41, 1);
		p3.add(lab4);
		lO3.setDimensions(0.43, 1);
		p3.add(average);
	}

	/**
	 * Sets all Components either to visible or unvisible.
	 * 
	 * @param b
	 *            boolean : true if visible, false otherwise.
	 */
	public void setCompVisible(boolean b) {
		subCat.setVisible(b);
		gammaBox.setVisible(b);
		alpha.setVisible(b);
		average.setVisible(b);
		lab1.setVisible(b);
		lab2.setVisible(b);
		lab3.setVisible(b);
		lab4.setVisible(b);
	}

	public void actionPerformed(ActionEvent e) {
		if (gammaBox.getSelectedItem().toString().equals("fixed")) {
			alpha.setEnabled(true);
			alpha.requestFocus();
		} else {
			alpha.setEnabled(false);
		}
	}

	class CustomTextField extends JTextField {
		/**
		 * default id
		 */
		private static final long serialVersionUID = 1L;
		private boolean type;

		/**
		 * Constructor method.
		 * 
		 * @param string
		 *            String : double value to set as default in text field.
		 * @param type
		 *            boolean : whether it is an integer or double. True if
		 *            Double value, false if not.
		 */
		public CustomTextField(String string, boolean type) {
			super(string);
			this.type = type;
			if (type) {
				if (!isDouble(string)) {
					this.setText("");
				}
			} else {
				if (!isInteger(string)) {
					this.setText("");
				}
			}
		}

		/**
		 * Tests whether the input String is convertable into an Integer.
		 * 
		 * @param str
		 *            String : input String.
		 * @return 
		 * boolean : true if convertible, false otherwise.
		 */
		private boolean isInteger(String str) {
			try {
				Integer.parseInt(str);
				return true;
			} catch (NumberFormatException e) {
				return false;
			}
		}

		@Override
		protected Document createDefaultModel() {
			return new doubleOnlyDocument();
		}

		/**
		 * Tests whether the input String is convertable into an Double.
		 * 
		 * @param str
		 *            String : input String.
		 * @return 
		 * boolean : true if convertible, false otherwise.
		 */
		private boolean isDouble(String str) {
			try {
				Double.parseDouble(str);
				return true;
			} catch (NumberFormatException e) {
				if (str.equals(".") && !getGammaString().contains(".")) {
					return true;
				}
				return false;
			}
		}

		/**
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
				if (type) {
					if (str == null || !isDouble(str)) {
						return;
					}
					super.insertString(offs, str, a);
				} else {
					if (str == null || !isInteger(str)) {
						return;
					}
					super.insertString(offs, str, a);
				}
			}
		}
	}

	/**
	 * Retrieves the String of the Gamma value from the Costume TextField.
	 * 
	 * @return 
	 * String : Gama value as a String.
	 */
	public String getGammaString() {
		return alpha.getText();
	}

	/**
	 * Retrieves the number of subCategories as a String.
	 * 
	 * @return 
	 * String : "4" by default otherwise what the user typed in.
	 */
	public String getNumSubCat() {
		if (subCat.getText().equals("")) {
			return 4 + "";
		}
		return subCat.getText();
	}

	/**
	 * Retrieves the Alpha value set by the user.
	 * 
	 * @return 
	 * String : returns "e" for estimated and the Alpha value as
	 * specified by the user otherwise.
	 */
	public String getAlpha() {
		if (gammaBox.getSelectedItem().toString().equals("estimated")) {
			return "e";
		} else if (!alpha.getText().equals("")) {
			return alpha.getText();
		}
		return "e";
	}

	/**
	 * Retrieves what average type is selected.
	 * 
	 * @return 
	 * String : "mean" or "median".
	 */
	public String getUseMedian() {
		return average.getSelectedItem().toString();
	}

	public void focusGained(FocusEvent e) {
	}

	public void focusLost(FocusEvent e) {
		if (alpha.getText().equals("")) {
			gammaBox.setSelectedIndex(0);
		}
	}
}
